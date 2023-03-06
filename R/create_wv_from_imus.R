library(imudata)
library(wv)
data(imu6)
MTiG = list(
  Gyro.X = wvar(imu6[,1]),
  Gyro.Y = wvar(imu6[,2]),
  Gyro.Z = wvar(imu6[,3]),
  Accel.X = wvar(imu6[,4]),
  Accel.Y = wvar(imu6[,5]),
  Accel.Z = wvar(imu6[,6]),
  freq = 100)

data(navchip)

navchip = list(
  Gyro.X = wvar(navchip[,1]),
  Gyro.Y = wvar(navchip[,2]),
  Gyro.Z = wvar(navchip[,3]),
  Accel.X = wvar(navchip[,4]),
  Accel.Y = wvar(navchip[,5]),
  Accel.Z = wvar(navchip[,6]),
  freq = attributes(navchip)$freq)

data(imar.gyro)

imar = list(
  Gyro.X = wvar(imar.gyro[,1]),
  Gyro.Y = wvar(imar.gyro[,2]),
  Gyro.Z = wvar(imar.gyro[,3]),
  freq = attributes(imar.gyro)$freq)

data("ln200.gyro")
data("ln200.accel")

ln200 = list(
  Gyro.X = wvar(ln200.gyro[,1]),
  Gyro.Y = wvar(ln200.gyro[,2]),
  Gyro.Z = wvar(ln200.gyro[,3]),
  Accel.X = wvar(ln200.accel[,1]),
  Accel.Y = wvar(ln200.accel[,2]),
  Accel.Z = wvar(ln200.accel[,3]),
  freq = attributes(ln200.gyro)$freq)

data = list(MTiG = MTiG,
            navchip = navchip,
            imar = imar,
            ln200 = ln200)


library(gmwm)
fit = gmwm(WN()+ RW(), input = data[[3]][[1]])
fit
str(data[[2]][[1]])
save(data, file = "data/imudata.RData")
