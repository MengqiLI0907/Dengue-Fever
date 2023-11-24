library(dlnm) ; library(mgcv) ; library(splines) ;
library(tsModel); library(reshape);library(ggplot2);library(dplyr);
library(mice)


data <- read.csv("F:\\GeoHealth_NUS\\Data\\State2017-2022\\State_Monthly2017-2022.csv")
# to analysis specific part
# 使用is.na()函数找到缺失值的位置
# missing_indices <- which(is.na(data$cases))
# print(missing_indices)
# # 后向填充（使用后一个非缺失值进行填补）
# for (i in rev(missing_indices)) {
#   data$cases[i] <- data$cases[i + 1]
# }
data[data$cases>9999,'cases'] <- NA  # 这里比较奇怪，pm10中有缺失值所以不能赋值为NA，插补后就可以了，没太弄懂哈哈哈
miceMod <- mice(data)  # 进行mice插值
data <- complete(miceMod)  # 生成完整数据
# data <- filter(data,Country=='MY' | Country=='SG')
data <- filter(data,Country=='MY' | Country=='SG')%>%filter(year<=2019)
# data <- filter(data,Country=='TH')%>%filter(year>2019)
# data <- filter(data, year>2019)
# set the lag and free-degree
lag_parameter <- 3
degree_parameter <- 2

# 降水量分析

# precipitation median value
knot.pr = round(median(data$pr))
# 本次建模使用了特征维度使用了ns，滞后维度使用了poly
cb.pr = crossbasis(data$pr, lag=lag_parameter, argvar=list(fun="ns", knots= knot.pr),
                   arglag=list(fun="poly",degree=degree_parameter))
# 拟合类 poission 模型
model1 = glm(cases ~ cb.pr + ns(time,5*1)+month,family=quasipoisson('log'), data)

# 拟合模型计算 RR 以及 RR 的置信区间， cen 选取对照点作为 RR 分母，bylag 表示 lag 的切片单位。
pred1.pr = crosspred(cb.pr, model1, cen=knot.pr, bylag=0.2, cumul=TRUE)

# 暴露-反应关系
plot(pred1.pr, "contour", xlab="precipitation", key.title=title("RR"),
     plot.title=title(xlab="precipitation (mm)",ylab="Lag (months)",cex.main=2,cex.lab=1.5))

plot(pred1.pr,"slices", lag=list(1,2,3),xlab="precipitation(mm)", ylab="RR",col=3,ci.arg=list(density=15,lwd=2))

crall <- crossreduce(cb.pr,model1, cen=knot.pr,type="overall",lag=c(1,3))
plot(crall,xlab="precipitation(mm)",ylab="RR",lwd=1.5)

mtext(text="Overall cumulative association 1-3 months",cex=0.89)

# 相对湿度分析
knot.rh = round(median(data$rh))

cb.rh = crossbasis(data$rh, lag=lag_parameter, argvar=list(fun="ns", knots= knot.rh),
                   arglag=list(fun="poly",degree=degree_parameter)) # 本次建模使用了特征维度使用了ns，滞后维度使用了poly
model2 = glm(cases ~ cb.rh + ns(time,5*1)+month,family=quasipoisson(), data) # 拟合类 poission 模型

pred1.rh = crosspred(cb.rh, model2, cen=knot.rh, bylag=0.2) # 拟合模型计算 RR 以及 RR 的置信区间， cen 选取对照点作为 RR 分母，bylag 表示 lag 的切片单位。

plot(pred1.rh, "contour", xlab="relative humidity(%)", key.title=title("RR"),
     plot.title=title(xlab="relative humidity(%)",ylab="Lag (months)",cex.main=2,cex.lab=1.5))

plot(pred1.rh,"slices", lag=list(1,2,3), col=3, ylab="RR",xlab="relative humidity(%)", ci.arg=list(density=15,lwd=2))

crall <- crossreduce(cb.rh,model2,cen=knot.rh,type="overall",lag=c(1,3))
plot(crall,xlab="relative humidity(%)",ylab="RR",lwd=1.5)
mtext(text="Overall cumulative association 1-3 months",cex=0.89)


# 最大温度分析
knot.tmmx = round(median(data$tmmx))

cb.tmmx = crossbasis(data$tmmx, lag=lag_parameter, argvar=list(fun="ns", knots= knot.tmmx),
                     arglag=list(fun="poly",degree=degree_parameter)) # 本次建模使用了特征维度使用了ns，滞后维度使用了poly
model3 = glm(cases ~ cb.tmmx + ns(time,5*1)+month,family=quasipoisson(), data) # 拟合类 poission 模型

pred1.tmmx = crosspred(cb.tmmx, model3, cen=knot.tmmx, bylag=0.2) # 拟合模型计算 RR 以及 RR 的置信区间， cen 选取对照点作为 RR 分母，bylag 表示 lag 的切片单位。

plot(pred1.tmmx, "contour", xlab="TempMax(°C)", key.title=title("RR"),
     plot.title=title(xlab="TempMax",ylab="Lag (months)",cex.main=2,cex.lab=1.5))

plot(pred1.tmmx,"slices", lag=list(1,2,3), col=3, ylab="RR", xlab="TempMax(°C)",ci.arg=list(density=15,lwd=2))

crall <- crossreduce(cb.tmmx,model3,cen=knot.tmmx,type="overall",lag=c(1,3))
plot(crall,xlab="TempMax",ylab="RR")
mtext(text="Overall cumulative association 1-3 months",cex=0.89)


# 最低温度分析
knot.tmmn = round(median(data$tmmn))
cb.tmmn = crossbasis(data$tmmn, lag=lag_parameter, argvar=list(fun="ns", knots= knot.tmmn),
                     arglag=list(fun="poly",degree=degree_parameter)) # 本次建模使用了特征维度使用了ns，滞后维度使用了poly

model4 = glm(cases ~ cb.tmmn + ns(time,5*1)+month,family=quasipoisson('log'), data) # 拟合类 poission 模型

pred1.tmmn = crosspred(cb.tmmn, model4, cen=knot.tmmn, bylag=0.2) # 拟合模型计算 RR 以及 RR 的置信区间， cen 选取对照点作为 RR 分母，bylag 表示 lag 的切片单位。

plot(pred1.tmmn, "contour", xlab="TempMin", key.title=title("RR"),
     plot.title=title(xlab="TempMin",ylab="Lag (months)",cex.main=2,cex.lab=1.5))

plot(pred1.tmmn,"slices", lag=list(1,2,3), col=3, ylab="RR",xlab="TempMin(°C)", ci.arg=list(density=15,lwd=2))

crall <- crossreduce(cb.tmmn,model4,cen=knot.tmmn,type="overall",lag=c(1,3))
plot(crall,xlab="TempMin",ylab="RR",lwd=1.5)
mtext(text="Overall cumulative association 1-3 months")