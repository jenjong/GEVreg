if(!require(splines)){ install.packages('splines')}; require(splines)
load("C:\\Users\\UOS\\Dropbox\\Extreme value\\Pr_48.RData")
Pr_48 <- Pr_48[Pr_48$stnlds!=115,]  # stnlds : 지역번호
head(Pr_48)
x <- Pr_48$pr       # 강수량   
latbs <- bs(Pr_48$lat,df=4) # 위도
longbs <- bs(Pr_48$long,df=4) # 경도 
# 2-dimensional tensor product basis
tensorbs <- do.call('cbind', lapply(1:ncol(latbs), function(i) latbs[, i] * longbs)) 
