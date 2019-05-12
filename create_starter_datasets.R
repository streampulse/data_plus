library(dplyr)

sp = read.csv('~/Downloads/all_sp_data.csv', stringsAsFactors=FALSE)
flagsites = c('751MUD','BEC','BRW','BUNN','Cameo','ColeMill','Eno','Icacos',
    'ICHE2700','LV','Mud','MV','NeuseDown1','NeuseDown2A','NeuseDown2B',
    'NeuseImp1','NeuseImp2','NeuseUp1','NeuseUp2','NHC','NR1000','OC',
    'Potash','QS','SDW','SF2500','SF2800','SF700','Simms','Sisimiut','Stony',
    'UEno','UNHC','WB','WS1500','YRN42')

sp = filter(sp, siteID %in% flagsites)

write.csv(sp, '~/git/data+/starter_kit/flagged_sites.csv', row.names=FALSE)
