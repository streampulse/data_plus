library(dplyr)
library(zoo)

sp = read.csv('~/git/_data+/starter_kit/flagged_sites.csv',
    stringsAsFactors=FALSE)

#prepare data; isolate one sensor series from one site
foc_site = 'NHC'
subset = filter(sp, siteID == foc_site)
rm(sp); gc() #free up memory
subset = subset[! duplicated(subset),]
subset$dateTimeUTC = as.POSIXct(subset$dateTimeUTC, tz='UTC')

variables = unique(subset$variable)
variables = variables[variables != 'Battery_V']
foc_var = variables[5]
series = filter(subset, variable == foc_var) %>%
    arrange(dateTimeUTC)

#remove duplicate rows, preferring "bad data" flags (these will usually be outliers)
is_dupe = duplicated(series$dateTimeUTC) |
    duplicated(series$dateTimeUTC, fromLast=TRUE)
dupes = series[is_dupe,]
dupe_key = factor(paste(dupes$dateTimeUTC, dupes$value))
keepers = data.frame()
for(k in levels(dupe_key)){
    dupe_group_inds = which(dupe_key == k)
    dupe_group = dupes$flagID[dupe_group_inds]
    keep = dupe_group_inds[dupe_group == 'Bad Data'][1]
    if(is.na(keep)) keep = dupe_group_inds[1]
    keepers = rbind(keepers, dupes[keep,])
}

series = filter(series, ! is_dupe) %>%
    bind_rows(keepers) %>%
    arrange(dateTimeUTC)

#fill out missing sample times with NAs; visualize
sample_interval = as.numeric(difftime(series$dateTimeUTC[2],
    series$dateTimeUTC[1]))

full_series = seq.POSIXt(series$dateTimeUTC[1],
    series$dateTimeUTC[nrow(series)], by=paste(sample_interval, 'min'))
series = full_join(data.frame(dateTimeUTC=full_series), series)

plot(series$dateTimeUTC, series$value, type='l', ylab=foc_var)

#convert to time series object
minutes_per_day = 1440
samp_per_day = minutes_per_day / sample_interval
tm = series$value
tm = ts(tm, deltat = 1/samp_per_day)

#heuristic outlier detection method ####

#linearly interpolate NAs
tm = na.approx(tm)
par(mfrow=c(3,1))
plot(series$dateTimeUTC, tm, col='orange', ylab=foc_var, type='l', xlab='',
    bty='l')
lines(series$dateTimeUTC, series$value, col='gray')
legend('bottomleft', legend='interpolated', lty=1, col='orange', bty='n')

#locate gobal outliers (>4sd beyond the mean)
outlier_list = list()

ts_u = mean(tm, na.rm=TRUE)
ts_sd = sd(tm, na.rm=TRUE)
big_outliers_h = which(tm > ts_u + (4 * ts_sd))
big_outliers_l = which(tm < ts_u - (4 * ts_sd))
big_outliers = unique(c(big_outliers_h, big_outliers_l))

points(series$dateTimeUTC[big_outliers], tm[big_outliers], col='red', pch=20)
legend('bottom', legend='+/- 4 SD', pch=20, col='red', bty='n')

#get diffs in value between each time point and the next
diffs = diff(tm)

#get mean and sd of diffs (large ones are later referred to as "jumps")
#ignore diffs made by big outliers when calculating these statistics.
bo_diffs = union(big_outliers, big_outliers - 1)
bo_diffs = bo_diffs[bo_diffs > 0]
u = mean(diffs[-bo_diffs], na.rm=TRUE)
sd = sd(diffs[-bo_diffs], na.rm=TRUE)

#expand sd threshold until <3% of diffs are outside sd bounds
#these are now potential outliers in addition to the globals above
sd_scaler = 1.8
big_jump_prop = Inf
if(! is.na(sd)){
    while(big_jump_prop > 0.03){
        sd_scaler = sd_scaler + 0.2
        big_jump_prop = sum(diffs > sd_scaler * sd) / length(diffs)
    }
}

#get indices of large jumps between successive points
pos_jumps = which(diffs > sd_scaler * sd)
neg_jumps = which(diffs < -sd_scaler * sd)
jump_inds = sort(c(pos_jumps, neg_jumps))

plot(diffs, col='gray', ylab=foc_var, main='diffs', xlab='', bty='l')
# xlim=c(as.POSIXct('2017-04-05 18:00:00'),
#     as.POSIXct('2017-04-05 20:00:00')))
# xlim=c(203.2,203.4),
# ylim=c(-50, 50))
points(time(tm)[pos_jumps + 1], diffs[pos_jumps], col='orange1')
points(time(tm)[neg_jumps + 1], diffs[neg_jumps], col='orange3')
legend('topleft', legend=c('positive jumps', 'negative jumps'),
    pch=1, col=c('orange1', 'orange3'), bty='n')

#if there are no such jumps to speak of, grab global outliers
#(or nothing) and move on
if(length(pos_jumps) == 0 | length(neg_jumps) == 0){
    if(length(big_outliers)){
        outlier_list[[col]] = big_outliers
    } else {
        outlier_list[[col]] = 'NONE'
    }
    names(outlier_list)[col] = colnames(df)[col]
    stop('done')
}

#an outlier as defined here must consist of a pair of positive and
#negative jumps. here we get the run lengths of consecutive positive
#(1) and negative (0) jumps
runs = rle(as.numeric(jump_inds %in% pos_jumps))
ends = cumsum(runs$lengths)
runs = cbind(values=runs$values, starts=c(1, ends[-length(ends)] + 1),
    stops=ends, lengths=runs$lengths, deparse.level=1)

lr = runs[,'lengths'] > 3 #runs > 3 are considered "long"
long_runs = runs[lr, 2:3, drop=FALSE]

#and then decide which ones to "keep", i.e. which ones may be actual outliers
keep = numeric()
if(length(long_runs)){
    for(i in 1:nrow(long_runs)){
        r = long_runs[i,]

        #get the indices of the first and second jumps in the run,
        #and also the last and second-to-last
        j = jump_inds[unique(c(r[1], r[1] + 1, r[2] - 1, r[2]))]
        # points(time(tm)[j], diffs[j], col='blue')

        #if both pairs' values are quite close to each other in time,
        #assume this is one long run of related jumps, e.g. the rising
        #limb of a storm, i.e. a true data feature and not an outlier
        if(abs(j[2]-j[1]) < 8 & abs(j[length(j)]-j[length(j)-1]) < 8){
            next
        }

        #otherwise it's a run of real jumps (a real data
        #feature) followed by a potential outlier, or vice-versa

        #find out which end of the run has the shorter time interval
        #between adjacent jumps and assume that end does not contain
        #the potential outlier (crude)
        t = time(tm)[j + 1] #+1 to convert from diff inds to ts inds
        left_jump_interv = t[2] - t[1]
        right_jump_interv = t[length(t)] - t[length(t)-1]
        not_outlier = which.min(c(left_jump_interv, right_jump_interv))

        #store the index of the other end of the run. this one might be
        #an actual outlier (a positive-negative or negative-positive jump pair,
        #hereinafter referred to as a "+/- pair" or posNeg_jump_pair)
        keep = append(keep, i[-not_outlier])
        # keep = append(keep, ifelse(not_outlier == 2, j[1], j[length(j)]))
    }
}

#all remaining inds may represent +/- pairs
short_run_inds = as.vector(runs[! lr, 2:3])
posNeg_jump_pairs = sort(unique(c(keep, short_run_inds)))
outlier_inds = jump_inds[posNeg_jump_pairs]

# #this section just for testing
# regulars = jump_inds[sort(unique(short_run_inds))]
# keepers = jump_inds[sort(unique(keep))]
# lra = as.vector(t(runs[lr, 2:3]))
# gfgf = mapply(seq, lra[seq(1, length(lra), 2)], lra[seq(2, length(lra), 2)])
# rbc = rainbow(10)
# rbcr = rep(rbc, length.out=length(gfgf))
# cols = rep(rbcr, times=sapply(gfgf, length))
# alllr = jump_inds[sort(unlist(gfgf))]
#
#
# plot(diffs, col='gray', ylab=foc_var, main='diffs', xlab='', bty='l', type='l',
#     # xlim=c(as.POSIXct('2017-04-05 18:00:00'),
#     #     as.POSIXct('2017-04-05 20:00:00')))
#     xlim=c(203.2,203.4),
#     # xlim=c(395.4,395.5),
#     # xlim=c(405.2,405.3),
#     # xlim=c(195,205),
#     # xlim=c(213,216),
#     # xlim=c(1.5,2),
#     # ylim=c(-5, 5))
#     # ylim=c(-.5, .5))
#     ylim=c(-50, 50))
# points(time(tm)[outlier_inds + 1], diffs[outlier_inds], col='green')
# points(time(tm)[alllr + 1], diffs[alllr], col=cols)
# points(time(tm)[keepers + 1], diffs[keepers], col='black')
# points(time(tm)[regulars + 1], diffs[regulars], col='orange3')


#winnow them down using various heuristics
n_outlier_pieces = Inf #an outlier piece is one unidirectional jump
counter = 0 #dont loop for too long
removed_first = removed_last = FALSE
if(length(outlier_inds) == 1){
    outlier_ts = 'NONE' #if only one ind, can't be a +/- pair
} else {
    while(length(outlier_inds) > 1 & n_outlier_pieces > 50 & counter < 6){

        outdif = diff(outlier_inds)
        rm_multjump = rm_oneway = NULL

        #now "jump" refers to the gap between potential outlier indices
        short_jumps = outdif < 15

        #***if the first jump is long, assume the first outlier piece
        #is a component of a data feature that got cut off by the border
        #of the time window, (which will cause a "frameshift mutation"
        #downstream, so remove it).
        #This step is only relevant when this script is used in production.
        #for Data+ we're operating on a full time series, not a moving window,
        #so it shouldn't be necessary.
        if(! short_jumps[1]){ # && ! removed_first){
            outlier_inds = outlier_inds[-1]
            removed_first = TRUE
            # counter = counter + 1
            next
        }


        #this works just like the *** comment above, but on the tail
        #of the series. sorry for the lack of naming consistency here.
        length(outlier_inds)
        big_outdif = outdif > 15
        if(big_outdif[length(big_outdif)]){ # && ! removed_last){
            outlier_inds = outlier_inds[-length(outlier_inds)]
            removed_last = TRUE
            # counter = counter + 1
            next
        }

        #a multijump is a series of jumps in the same direction within
        #a short duration. not likely to be a real outlier.
        jump_runs = rle(as.numeric(short_jumps))
        ends2 = cumsum(jump_runs$lengths)
        jump_runs = cbind(values=jump_runs$values,
            starts=c(1, ends2[-length(ends2)] + 1),
            stops=ends2, lengths=jump_runs$lengths, deparse.level=1)

        multijumps = jump_runs[jump_runs[,'lengths'] > 1, 2:3, drop=FALSE]
        if(length(multijumps)){

            #filter jumps with long gaps between (keep only short multijumps)
            multijumps = multijumps[short_jumps[multijumps[,'starts']],,
                drop=FALSE]

            #interpolate indices within multijumps that do not
            #themselves represent jumps
            seq_list = mapply(function(x, y){ seq(x, y + 1, 1) },
                multijumps[,1], multijumps[,2],
                SIMPLIFY=FALSE)

            #these are the indices of mutijumps, which will be ignored as
            #potential outliers
            rm_multjump = unlist(seq_list)
        }

        #if none of the outlier pieces is near another, there can be
        #no +/- pairs, so assume no outliers
        if(all(big_outdif)){
            outlier_ts = 'NONE'
            break
        } else {

            #find remaining one-way jumps (outlier pieces)
            same_sign_runs = rle(as.numeric(big_outdif))
            same_sign_run_ends = cumsum(same_sign_runs$lengths)
            same_sign_run_starts = c(1,
                same_sign_run_ends[-length(same_sign_run_ends)] + 1)

            if(length(same_sign_runs)){
                l = same_sign_runs$lengths
                one_way_jumps = numeric()
                for(i in 1:length(l)){
                    if(l[i] > 1){
                        s = same_sign_run_starts[i]
                        one_way_jumps = append(one_way_jumps,
                            seq(s, s + (l[i] - 2), 1))
                    }
                }
            }

            #just a bit of hacky filtering to isolate the one-way jumps
            if(length(one_way_jumps)){
                one_way_jumps = one_way_jumps[big_outdif[one_way_jumps]]
                outdif_filt = outdif
                outdif_filt[-one_way_jumps] = 0
                rm_oneway = which(outdif_filt > 0) + 1
            }
        }
        removals = unique(c(rm_multjump, rm_oneway))
        if(!is.null(removals)){
            outlier_inds = outlier_inds[-removals]
        }

        #if there's an even number of outlier pieces remaining,
        #interpolate between the members of each pos/neg pair and
        #contribute to list of outlier components to be identified
        #for flagging in qa/qc
        if(length(outlier_inds) %% 2 == 0){
            outlier_inds = matrix(outlier_inds, ncol=2, byrow=TRUE)
            seq_list = mapply(function(x, y){ seq(x+1, y, 1) },
                outlier_inds[,1], outlier_inds[,2],
                SIMPLIFY=FALSE)
            outlier_ts = unlist(seq_list)

            #otherwise get rid of the least extreme outlier piece
        } else {
            smallest_diff = which.min(abs(diffs[outlier_inds]))
            outlier_inds = outlier_inds[-smallest_diff]
            counter = counter + 1
            next
        }
        # print(length(outlier_inds))
        # print(counter)
        n_outlier_pieces = length(outlier_ts)
        counter = counter + 1
    }
}

#if the while loop completes and there are still tons of pieces,
#assume there was a frameshift or something and abort.
#the bit about outlier_ts existing is for compatibility once this R
#code gets translated to python (which can't handle R's NULL)
if(n_outlier_pieces > 50 || !exists('outlier_ts') ||
        is.null(outlier_ts)){ #deal with R-py null value mismatch
    outlier_ts = 'NONE'
}

#bring in the global outliers from above
if(length(outlier_ts) == 1 && outlier_ts == 'NONE'){
    if(length(big_outliers)){
        outlier_ts = unique(big_outliers)
    } else {
        outlier_ts = 'NONE'
    }
} else {
    outlier_ts = unique(c(outlier_ts, big_outliers))
}

#list appending is a relic from production, where this runs in a loop
outlier_list[[1]] = outlier_ts
names(outlier_list)[1] = colnames(df)[1]

#visualize global and local potential outliers together
plot(series$dateTimeUTC, tm, col='orange', ylab=foc_var, type='l', xlab='',
    bty='l')
    # xlim=c(as.POSIXct('2017-04-05 18:00:00'),
    #     as.POSIXct('2017-04-05 20:00:00')),
    # # bty='l', xlim=c(as.POSIXct('2017-04-05'), as.POSIXct('2017-04-07')),
    # ylim=c(0, 50))
lines(series$dateTimeUTC, series$value, col='gray')
points(series$dateTimeUTC[outlier_ts], tm[outlier_ts], col='darkred', pch=20)
points(series$dateTimeUTC[big_outliers], tm[big_outliers], col='red', pch=20)
legend('bottom', legend=c('global outlier', 'local outlier'),
    pch=20, col=c('red', 'darkred'), bty='n')
