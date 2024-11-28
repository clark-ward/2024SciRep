# Associated with Ward et al. (2024/5) "Enamel carbon, oxygen, and strontium isotopes reveal limited mobility in an extinct rhinoceros at Ashfall Fossil Beds, Nebraska, USA" _Scientific Reports_

dir <- c("/Users/cward747/Desktop/R/2024SciRep")
setwd(dir)

# install packages
{
  install.packages("PMCMRplus")
}

# load packages
{
  library("PMCMRplus")
}

# load in data 
{
  bulk_data <- read.csv("data/Bulk_data.csv", header = T, sep = ",", stringsAsFactors = F, fileEncoding = "UTF-8")
  serial_data <- read.csv("data/Serial_data.csv", header = T, sep = ",", stringsAsFactors = F, fileEncoding = "UTF-8")
  AFB_data <- read.csv("data/AFB_data.csv", header = T, sep = ",", stringsAsFactors = F, fileEncoding = "UTF-8")
}

# data manipulation
{
  # separating things out for brevity later
  male.bulk <- bulk_data[bulk_data$sex == "m",]
  female.bulk <- bulk_data[bulk_data$sex == "f",]
  m2.bulk <- bulk_data[bulk_data$tooth.position == "m2",]
  m3.bulk <- bulk_data[bulk_data$tooth.position == "m3",]
  m2.serial <- serial_data[serial_data$tooth.position == "m2",]
  m3.serial <- serial_data[serial_data$tooth.position == "m3",]
  Woofy.results <- serial_data[serial_data$lab.ID == "AFB-02",]
  John.results <- serial_data[serial_data$lab.ID == "AFB-25",]
  Mary.results <- serial_data[serial_data$lab.ID == "AFB-07",]
  Sparky.results <- serial_data[serial_data$lab.ID == "AFB-30",]
  Longi <- AFB_data[AFB_data$family == "Blastomerycinae",]
  Horses <- AFB_data[AFB_data$family == "Equidae",]
  Camels <- AFB_data[AFB_data$family == "Camelidae",]
  
  #calculating m3-2 differences, or capital delta (cap del)
  cap.del <- m2.bulk
  cap.del$DSr.iso <- m3.bulk$Sr.iso - m2.bulk$Sr.iso
  cap.del$Dd13C <- m3.bulk$d13C - m2.bulk$d13C
  cap.del$Dd18O.VPDB <- m3.bulk$d18O.VPDB - m2.bulk$d18O.VPDB
  cap.del$Dd18O.VSMOW <- m3.bulk$d18O.VSMOW - m2.bulk$d18O.VSMOW
  cap.del$grp <- c() ; cap.del[cap.del$sex == "m", "grp"] <- 1 ; cap.del[cap.del$sex == "f", "grp"] <- 2
  
  #calculating m2+m3 average for family-level comparisons
  AFB_data_2 <- AFB_data[AFB_data$taxon != "Teleoceras major",]
  AFB_data_2$sample.no <- NA
  AFB_data_2$grp <- NA
  avg_bulk <- m2.bulk
  avg_bulk$family <- "Rhinocerotidae"
  avg_bulk$tooth.position <- "average"
  avg_bulk$color <- NA
  avg_bulk$pch <- NA
  avg_bulk$Sr.iso <- (m2.bulk$Sr.iso + m3.bulk$Sr.iso) / 2
  avg_bulk$d13C <- (m2.bulk$d13C + m3.bulk$d13C) / 2
  avg_bulk$d18O.VPDB <- (m2.bulk$d18O.VPDB + m3.bulk$d18O.VPDB) / 2
  avg_bulk$d18O.VSMOW <- (m2.bulk$d18O.VSMOW + m3.bulk$d18O.VSMOW) / 2
  AFB_data_2 <- rbind(AFB_data_2, avg_bulk)
}

# objects to make plotting easier
{
  VPDB.lim <- c(-10, -0.8)
  color.list <- c("darkorange", "deepskyblue", "indianred3", "aquamarine3", "paleturquoise4", "firebrick4", "lightgoldenrod4", "khaki3", "lightskyblue3", "olivedrab", "saddlebrown", "lavender", "slateblue1")
  Sr.iso.uncertainty <- 0.00003
  C.lab <- expression(bold(delta^13*"C"["VPDB"]* " (\u2030)"))
  O.lab <- expression(bold(delta^18*"O"["VPDB"]* " (\u2030)"))
  Sr.lab <- expression(bold(""^87*"Sr/"^86*"Sr"))
}

# non-parametric analyses
{
# m2 vs m3 PAIRED
  # males
  wilcox.test(male.bulk[male.bulk$tooth.position == "m2", "d13C"], male.bulk[male.bulk$tooth.position == "m3", "d13C"], paired = T)
  fligner.test(male.bulk$d13C, male.bulk$tooth.position)
  wilcox.test(male.bulk[male.bulk$tooth.position == "m2", "d18O.VPDB"], male.bulk[male.bulk$tooth.position == "m3", "d18O.VPDB"], paired = T)
  fligner.test(male.bulk$d18O.VPDB, male.bulk$tooth.position)
  wilcox.test(male.bulk[male.bulk$tooth.position == "m2", "Sr.iso"], male.bulk[male.bulk$tooth.position == "m3", "Sr.iso"], paired = T)
  fligner.test(male.bulk$Sr.iso, male.bulk$tooth.position)
  
  # females
  wilcox.test(female.bulk[female.bulk$tooth.position == "m2", "d13C"], female.bulk[female.bulk$tooth.position == "m3", "d13C"], paired = T)
  fligner.test(female.bulk$d13C, female.bulk$tooth.position)
  wilcox.test(female.bulk[male.bulk$tooth.position == "m2", "d18O.VPDB"], female.bulk[female.bulk$tooth.position == "m3", "d18O.VPDB"], paired = T)
  fligner.test(female.bulk$d18O.VPDB, female.bulk$tooth.position)
  wilcox.test(female.bulk[male.bulk$tooth.position == "m2", "Sr.iso"], female.bulk[male.bulk$tooth.position == "m3", "Sr.iso"], paired = T)
  fligner.test(female.bulk$Sr.iso, female.bulk$tooth.position)
  
# males vs females UNPAIRED
  # m2
  wilcox.test(m2.bulk[m2.bulk$sex == "m", "d13C"], m2.bulk[m2.bulk$sex == "f", "d13C"], paired = F)
  fligner.test(m2.bulk$d13C, m2.bulk$sex)
  wilcox.test(m2.bulk[m2.bulk$sex == "m", "d18O.VPDB"], m2.bulk[m2.bulk$sex == "f", "d18O.VPDB"], paired = F)
  fligner.test(m2.bulk$d18O.VPDB, m2.bulk$sex)
  wilcox.test(m2.bulk[m2.bulk$sex == "m", "Sr.iso"], m2.bulk[m2.bulk$sex == "f", "Sr.iso"], paired = F)
  fligner.test(m2.bulk$Sr.iso, m2.bulk$sex)
  
  # m3
  wilcox.test(m3.bulk[m3.bulk$sex == "m", "d13C"], m3.bulk[m3.bulk$sex == "f", "d13C"], paired = F)
  fligner.test(m3.bulk$d13C, m3.bulk$sex)
  wilcox.test(m3.bulk[m3.bulk$sex == "m", "d18O.VPDB"], m3.bulk[m3.bulk$sex == "f", "d18O.VPDB"], paired = F)
  fligner.test(m3.bulk$d18O.VPDB, m3.bulk$sex)
  wilcox.test(m3.bulk[m3.bulk$sex == "m", "Sr.iso"], m3.bulk[m3.bulk$sex == "f", "Sr.iso"], paired = F)
  fligner.test(m3.bulk$Sr.iso, m3.bulk$sex)
  
# differences UNPAIRED
  # between sexes
  wilcox.test(cap.del[cap.del$sex == "m", "Dd13C"], cap.del[cap.del$sex == "f", "Dd13C"], paired = F)
  fligner.test(cap.del$Dd13C, cap.del$sex)
  wilcox.test(cap.del[cap.del$sex == "m", "Dd18O.VPDB"], cap.del[cap.del$sex == "f", "Dd18O.VPDB"], paired = F)
  fligner.test(cap.del$Dd18O.VPDB, cap.del$sex)
  wilcox.test(cap.del[cap.del$sex == "m", "DSr.iso"], cap.del[cap.del$sex == "f", "DSr.iso"], paired = F)
  fligner.test(cap.del$DSr.iso, cap.del$sex)
  
# K-W, post hoc, and F-K tests BY FAMILY
  # carbon
  kruskal.test(formula = d13C ~ as.factor(family), data = AFB_data_2)
  summary(dscfAllPairsTest(formula = d13C ~ as.factor(family), data = AFB_data_2))
  fligner.test(formula = d13C ~ as.factor(family), data = AFB_data_2)
  
  # oxygen
  kruskal.test(formula = d18O.VPDB ~ as.factor(family), data = AFB_data_2)
  summary(dscfAllPairsTest(formula = d18O.VPDB ~ as.factor(family), data = AFB_data_2))
  fligner.test(formula = d18O.VPDB ~ as.factor(family), data = AFB_data_2)
  
  # strontium
  kruskal.test(formula = Sr.iso ~ as.factor(family), data = AFB_data_2)
  summary(dscfAllPairsTest(formula = Sr.iso ~ as.factor(family), data = AFB_data_2))
  fligner.test(formula = Sr.iso ~ as.factor(family), data = AFB_data_2)
}

# Figure 3
{
  cairo_pdf("figure 3.pdf", width = 9, height = 8)
  layout(matrix(c(1,2,3,4,5,6,7,7,7), ncol = 3, byrow = T), heights = c(3.5,3.5,1))
  par(mar = c(2,4,1.5,.5), mgp = c(2,1,0))
  
  #d13C
  boxplot(bulk_data$d13C ~ bulk_data$grp, ylab = expression(bold(delta^13*"C"["VPDB"]* " (\u2030)")), ylim = c(-10,-7.5), xaxt = "n", cex.lab = 1.5, whisklty = "solid", col = "white", cex.axis = 1.25, range = 0, xlab = NA)
  axis(side = 1, at = c(1:4), labels = as.vector(c("male", "male", "female", "female")), padj = -.6, cex.axis = 1.2)
  axis(side = 1, at = c(1:4), labels = as.vector(c(expression("M"[2]*"'s"), expression("M"[3]*"'s"),expression("M"[2]*"'s"),expression("M"[3]*"'s"))), padj = 1, cex.axis = 1.2)
  points(y = bulk_data[bulk_data$grp == 1, "d13C"], x = bulk_data$grp[1:5], col = color.list[1:5], cex = 2.5, pch = 15)
  points(y = bulk_data[bulk_data$grp == 2, "d13C"], x = bulk_data$grp[6:10], col = color.list[1:5], cex = 2.5, pch = 15)
  points(y = bulk_data[bulk_data$grp == 3, "d13C"], x = bulk_data$grp[11:18], col = c(color.list[6:11], "black", color.list[13]), bg = color.list[12], cex = 2.5, pch = c(16,16,16,16,16,16,21,16))
  points(y = bulk_data[bulk_data$grp == 4, "d13C"], x = c(4,4,3.95,3.95,4,4.05,4.05,4), col = c(color.list[6:11], "black", color.list[13]), bg = color.list[12], cex = 2.5, pch = c(16,16,16,16,16,16,21,16))
  segments(x0 = 1, x1 = 2, y0 = bulk_data[bulk_data$grp == 1, "d13C"], y1 = bulk_data[bulk_data$grp == 2, "d13C"], lty = 2, col = color.list[1:5])
  segments(x0 = 3, x1 = 4, y0 = bulk_data[bulk_data$grp == 3, "d13C"], y1 = bulk_data[bulk_data$grp == 4, "d13C"], lty = 2, col = color.list[6:13])
  
  #d18O
  boxplot(bulk_data$d18O.VPDB ~ bulk_data$grp, ylab = expression(bold(delta^18*"O"["VPDB"]* " (\u2030)")), ylim = c(-8,-4.5), xaxt = "n", cex.lab = 1.5, whisklty = "solid", col = "white", cex.axis = 1.25, xlab = NA)
  axis(side = 1, at = c(1:4), labels = as.vector(c("male", "male", "female", "female")), padj = -.6, cex.axis = 1.2)
  axis(side = 1, at = c(1:4), labels = as.vector(c(expression("M"[2]*"'s"), expression("M"[3]*"'s"),expression("M"[2]*"'s"),expression("M"[3]*"'s"))), padj = 1, cex.axis = 1.2)
  points(y = bulk_data[bulk_data$grp == 1, "d18O.VPDB"], x = c(.95,1,1,1.05,1), col = color.list[1:5], cex = 2.5, pch = 15)
  points(y = bulk_data[bulk_data$grp == 2, "d18O.VPDB"], x = bulk_data$grp[6:10], col = color.list[1:5], cex = 2.5, pch = 15)
  points(y = bulk_data[bulk_data$grp == 3, "d18O.VPDB"], x = c(3,3,2.95,3,3,3,3.05,3), col = c(color.list[6:11], "black", color.list[13]), bg = color.list[12], cex = 2.5, pch = c(16,16,16,16,16,16,21,16))
  points(y = bulk_data[bulk_data$grp == 4, "d18O.VPDB"], x = c(3.95,4,4.05,4,4,4,4,4), col = c(color.list[6:11], "black", color.list[13]), bg = color.list[12], cex = 2.5, pch = c(16,16,16,16,16,16,21,16))
  segments(x0 = 1, x1 = 2, y0 = bulk_data[bulk_data$grp == 1, "d18O.VPDB"], y1 = bulk_data[bulk_data$grp == 2, "d18O.VPDB"], lty = 2, col = color.list[1:5])
  segments(x0 = 3, x1 = 4, y0 = bulk_data[bulk_data$grp == 3, "d18O.VPDB"], y1 = bulk_data[bulk_data$grp == 4, "d18O.VPDB"], lty = 2, col = color.list[6:13])
  
  #Sr
  boxplot(bulk_data$Sr.iso ~ bulk_data$grp, ylim = c(0.70860,0.70875), ylab = expression(bold(""^87*"Sr/"^86*"Sr")), xaxt = "n", cex.lab = 1.5, whisklty = "solid", col = "white", cex.axis = 1.25, range = 0, xlab = NA)
  axis(side = 1, at = c(1:4), labels = as.vector(c("male", "male", "female", "female")), padj = -0.6, cex.axis = 1.2)
  axis(side = 1, at = c(1:4), labels = as.vector(c(expression("M"[2]*"'s"), expression("M"[3]*"'s"),expression("M"[2]*"'s"),expression("M"[3]*"'s"))), padj = 1, cex.axis = 1.2)
  points(y = bulk_data[bulk_data$grp == 1, "Sr.iso"], x = c(1,0.95,1,1,1.05), col = color.list[1:5], cex = 2.5, pch = 15)
  points(y = bulk_data[bulk_data$grp == 2, "Sr.iso"], x = bulk_data$grp[6:10], col = color.list[1:5], cex = 2.5, pch = 15)
  points(y = bulk_data[bulk_data$grp == 3, "Sr.iso"], x= c(rep(3, times = 6),2.95,3.05), col = c(color.list[6:11], "black", color.list[13]), bg = color.list[12], cex = 2.5, pch = c(16,16,16,16,16,16,21,16))
  points(y = bulk_data[bulk_data$grp == 4, "Sr.iso"], x = c(4,3.95,4.05,4,4,4,3.95,4.05), col = c(color.list[6:11], "black", color.list[13]), bg = color.list[12], cex = 2.5, pch = c(16,16,16,16,16,16,21,16))
  segments(x0 = 1, x1 = 2, y0 = bulk_data[bulk_data$grp == 1, "Sr.iso"], y1 = bulk_data[bulk_data$grp == 2, "Sr.iso"], lty = 2, col = color.list[1:5])
  segments(x0 = 3, x1 = 4, y0 = bulk_data[bulk_data$grp == 3, "Sr.iso"], y1 = bulk_data[bulk_data$grp == 4, "Sr.iso"], lty = 2, col = color.list[6:13])
  
  #Dd13C
  par(mar = c(2,4,2,.5), mgp = c(2,1,0))
  boxplot(cap.del$Dd13C ~ cap.del$grp, ylab = expression(bold(Delta*delta^13*"C"["M"["3"]*"-M"["2"]]* " (\u2030)")), ylim = c(-1, 1), xaxt = "n", cex.lab = 1.5, whisklty = "solid", col = "white", cex.axis = 1.25, xlab = NA)
  axis(side = 1, at = c(1:2), labels = as.vector(c("male", "female")), padj = -0.6, cex.axis = 1.2)
  points(y = cap.del[cap.del$sex == "m", "Dd13C"], x = rep(1, times = 5), col = color.list[1:5], cex = 2.5, pch = 15)
  points(y = cap.del[cap.del$sex == "f", "Dd13C"], x = c(1.95,2,2.05,1.95,2,1.9,2.1,2.05), col = c(color.list[6:11], "black", color.list[13]), bg = color.list[12], cex = 2.5, pch = c(16,16,16,16,16,16,21,16))
  
  #Dd18O
  boxplot(cap.del$Dd18O.VPDB ~ cap.del$grp, ylab = expression(bold(Delta*delta^18*"O"["M"["3"]*"-M"["2"]]* " (\u2030)")), ylim = c(-2, 2.5), xaxt = "n", cex.lab = 1.5, whisklty = "solid", col = "white", cex.axis = 1.25, xlab = NA)
  axis(side = 1, at = c(1:2), labels = as.vector(c("male", "female")), padj = -0.6, cex.axis = 1.2)
  points(y = cap.del[cap.del$sex == "m", "Dd18O.VPDB"], x = rep(1, times = 5), col = color.list[1:5], cex = 2.5, pch = 15)
  points(y = cap.del[cap.del$sex == "f", "Dd18O.VPDB"], x = c(rep(2, times = 6),1.95,2.05), col = c(color.list[6:11], "black", color.list[13]), bg = color.list[12], cex = 2.5, pch = c(16,16,16,16,16,16,21,16))
  
  #DSr
  boxplot(cap.del$DSr.iso ~ cap.del$grp, ylab = expression(bold(Delta*""^87*"Sr/"^86*"Sr"["M"["3"]*"-M"["2"]])), ylim = c(-6e-05, 6e-05), xaxt = "n", cex.lab = 1.5, whisklty = "solid", col = "white", cex.axis = 1.25, xlab = NA)
  axis(side = 1, at = c(1:2), labels = as.vector(c("male", "female")), padj = -0.6, cex.axis = 1.2)
  points(y = cap.del[cap.del$sex == "m", "DSr.iso"], x = rep(1, times = 5), col = color.list[1:5], cex = 2.5, pch = 15)
  points(y = cap.del[cap.del$sex == "f", "DSr.iso"], x = c(rep(2, times = 6),1.9,2.1), col = c(color.list[6:11], "black", color.list[13]), bg = color.list[12], cex = 2.5, pch = c(16,16,16,16,16,16,21,16))
  
  #legend
  par(mar = c(0,0.4,0,0.4))
  plot.new()
  legend("top", c("UNSM 52272 'Woofy'", "UNSM 27803 'Morris'", "UNSM 52236 'Bruce'", "UNSM 52288 'Fairfield'", "UNSM 27805 'Grandpa John'", "UNSM 27808 'Scarlet O'Hara'", "UNSM 52373 'Diller'", "UNSM 52286 'Mary Anning'", "UNSM 52282 'Jean'", "UNSM 27806 'Barb'", "UNSM 52289 'Joseph Leidy'", "UNSM 52269 'Pushy'", "UNSM 27807 'Sparky'", "Male (n=5)", "Female (n=8)"), col = c(color.list[1:11],"black",color.list[13],1,1), pt.bg = color.list[12], pch = c(15,15,15,15,15,16,16,16,16,16,16,21,16,15,16), pt.cex = 2.25, ncol = 3, text.width = 0.225, bty = "n", cex = 1.25)
  
  dev.off()
}

# This datum was excluded in summary statistics, so it is NA in the data file. I'm putting it in now for plotting. 
Mary.results[1,"Sr.iso"] <- 0.70856
# Mary.results[1,"Sr.iso"] <- NA; # Make sure it is NA for summary statistics!

# Figure 4
{
  # Function for each individual
  plot.serial <- function(data){
    plot(y= data$height, x= data$d18O.VPDB, xlim = VPDB.lim, col = "blue3", pch = 18, type = "o", ylim = c(1,45), ylab = NA, xlab = NA, cex = 1.25, yaxt = "n")
    points(y= data$height, x= data$d13C, col = "green4", pch = 24, type = "o", bg = "green4")
    par(new = T)
    plot(y= data$height, x= data$Sr.iso, ylim = c(1, 45), xlim = c(0.70835, 0.70875), xlab = NA, ylab = NA, axes = F, col = "darkorange", pch = 8, type = "o")
    mtext(expression( bold(delta^13*"C"* " & "* delta^18*"O"["VPDB"]* " (\u2030)")), side = 1, line = 2.6, cex = 1)
    axis(3, at = pretty(c(0.70855, 0.70875)), labels = T, tick = T)
    mtext(expression(bold(""^87*"Sr/"^86*"Sr")), side = 3, line = 2, cex = 1)
    legend("topleft", legend = c(expression(bold(delta^13*"C"["VPDB"]* " (\u2030)")), expression( bold(delta^18*"O"["VPDB"]* " (\u2030)")), expression(bold(""^87*"Sr/"^86*"Sr"))), cex = 1.2, col = c("green4", "blue3", "darkorange"), pch = c(17,18,8), pt.cex = c(1.5, 1.75, 1.5), bty = "n", horiz = T)
  }
  
  # Figure 4
  {
    cairo_pdf("figure 4.pdf", width = 7, height = 16)
    par(mfrow = c(4,2), cex.lab = 1, cex.main = 1.5)
    
    # Woofy's m2
    par(mar = c(4, 3.5, 4, 0.5))
    plot.serial(data = Woofy.results[Woofy.results$tooth.position == "m2",])
    axis(2, at = pretty(c(1,45)), labels = T, tick = T)
    mtext("Height from cervix (mm)", side = 2, line = 2, cex = 1)
    
    # Woofy's m3
    par(mar = c(4, 0.5, 4, 3.5))
    plot.serial(data = Woofy.results[Woofy.results$tooth.position == "m3",])
    axis(4, at = pretty(c(1,45)), labels = T, tick = T)
    mtext("Height from cervix (mm)", side = 4, line = 2, cex = 1)
    
    # Grandpa John's m2
    par(mar = c(4, 3.5, 4, 0.5))
    plot.serial(data = John.results[John.results$tooth.position == "m2",])
    axis(2, at = pretty(c(1,45)), labels = T, tick = T)
    mtext("Height from cervix (mm)", side = 2, line = 2, cex = 1)
    
    # Grandpa John's m3
    par(mar = c(4, 0.5, 4, 3.5))
    plot.serial(data = John.results[John.results$tooth.position == "m3",])
    axis(4, at = pretty(c(1,45)), labels = T, tick = T)
    mtext("Height from cervix (mm)", side = 4, line = 2, cex = 1)
    
    # Mary Anning's m2
    par(mar = c(4, 3.5, 4, 0.5))
    plot.serial(data = Mary.results[Mary.results$tooth.position == "m2",])
    axis(2, at = pretty(c(1,45)), labels = T, tick = T)
    mtext("Height from cervix (mm)", side = 2, line = 2, cex = 1)
    # Remember that the basal-most Sr value for this tooth was the anomalous datum and was excluded from analyses. I change how this point is plotted post-R. 
    
    # Mary Anning's m3
    par(mar = c(4, 0.5, 4, 3.5))
    plot.serial(data = Mary.results[Mary.results$tooth.position == "m3",])
    axis(4, at = pretty(c(1,45)), labels = T, tick = T)
    mtext("Height from cervix (mm)", side = 4, line = 2, cex = 1)
    
    # Sparky's m2
    par(mar = c(4, 3.5, 4, 0.5))
    plot.serial(data = Sparky.results[Sparky.results$tooth.position == "m2",])
    axis(2, at = pretty(c(1,45)), labels = T, tick = T)
    mtext("Height from cervix (mm)", side = 2, line = 2, cex = 1)
    
    # Sparky's m3
    par(mar = c(4, 0.5, 4, 3.5))
    plot.serial(data = Sparky.results[Sparky.results$tooth.position == "m3",])
    axis(4, at = pretty(c(1,45)), labels = T, tick = T)
    mtext("Height from cervix (mm)", side = 4, line = 2, cex = 1)
    
    dev.off()
  }
}

# Figure 5
{
   #yes this is just a big ugly block of code; don't judge me
    cairo_pdf("figure 5.pdf", width = 6, height = 13)
    par(mar = c(3, 5, 1, 1), cex = 1.2, cex.axis = 1.2, mfrow = c(3, 1))
    #d13C
    {
      boxplot(m2.serial$d13C, m3.serial$d13C, m2.bulk$d13C, m3.bulk$d13C, Longi$d13C, Camels$d13C, Horses$d13C, range = 0, ylim = c(-10.2, -6.7), xlab = NA, ylab = C.lab, col = "white", whisklty = "solid", xaxt = "n", cex.lab = 1.5)
      axis(side = 1, at = c(1:4), labels = as.vector(c(expression("M"[2]*"'s"), expression("M"[3]*"'s"),expression("M"[2]*"'s"),expression("M"[3]*"'s"))))
      points(Woofy.results$x, Woofy.results$d13C, pch = 22, bg = "green4", col = "white", cex = 2)
      points(John.results$x, John.results$d13C, pch = 22, bg = "white", col = "green4", cex = 2)
      points(Mary.results$x, Mary.results$d13C, pch = 21, bg = "green4", col = "white", cex = 2)
      points(Sparky.results$x, Sparky.results$d13C, pch = 21, bg = "white", col = "green4", cex = 2)
      points(male.bulk$x, male.bulk$d13C, pch = 22, bg = "green4", col = "black", cex = 2)
      points(female.bulk$x, female.bulk$d13C, pch = 21, bg = "green4", col = "black", cex = 2)
      points(c(5, 5, 5), Longi$d13C, pch = 9, col = "green4", cex = 2)
      points(c(6, 6), Camels$d13C, pch = 4, col = "green4", cex = 2)
      points(rep(7, length = 9), Horses$d13C, pch = 17, col = "green4", cex = 2)
      legend("topleft", legend = c(expression(bolditalic("Teleoceras major")), expression(bold("Serially sampled")), "UNSM 52272 'Woofy'", "UNSM 27805 'Grandpa John'", "UNSM 52286 'Mary Anning'", "UNSM 27807 'Sparky'", NA, expression(bold("Bulk sampled")),"Males (N = 5)", "Females (N = 8)", NA, NA), col = c(NA, NA, "white", "green4", "white", "green4", NA, NA, "black", "black", NA, NA), pt.bg = c(NA, NA, "green4", "white", "green4", "white", NA, NA, "green4", "green4", NA, NA), pch = c(NA, NA, 22, 22, 21, 21, NA, NA, 22, 21, NA, NA), cex = 1, ncol = 2, pt.cex = 2, text.width = c(2.25, 1.25))
      legend("bottomright", legend = c(expression(bold("Co-occurring ungulates")), "Blastomerycinae (N = 3)", "Camelidae (N = 2)", "Equidae (N = 9)"), col = c(NA, "green4", "green4", "green4"), pch = c(NA, 9, 4, 17), cex = 1, ncol = 1, pt.cex = 2)
    }
    #d18O.VPDB
    {
      boxplot(m2.serial$d18O.VPDB, m3.serial$d18O.VPDB, m2.bulk$d18O.VPDB, m3.bulk$d18O.VPDB, Longi$d18O.VPDB, Camels$d18O.VPDB, Horses$d18O.VPDB, range = 0, ylim = c(-8.1, -0.9), xlab = NA, ylab = O.lab, col = "white", whisklty = "solid", xaxt = "n", cex.lab = 1.5)
      axis(side = 1, at = c(1:4), labels = as.vector(c(expression("M"[2]*"'s"), expression("M"[3]*"'s"),expression("M"[2]*"'s"),expression("M"[3]*"'s"))))
      points(Woofy.results$x, Woofy.results$d18O.VPDB, pch = 22, bg = "blue3", col = "white", cex = 2)
      points(John.results$x, John.results$d18O.VPDB, pch = 22, bg = "white", col = "blue3", cex = 2)
      points(Mary.results$x, Mary.results$d18O.VPDB, pch = 21, bg = "blue3", col = "white", cex = 2)
      points(Sparky.results$x, Sparky.results$d18O.VPDB, pch = 21, bg = "white", col = "blue3", cex = 2)
      points(male.bulk$x, male.bulk$d18O.VPDB, pch = 22, bg = "blue3", col = "black", cex = 2)
      points(female.bulk$x, female.bulk$d18O.VPDB, pch = 21, bg = "blue3", col = "black", cex = 2)
      points(c(5, 5, 5), Longi$d18O.VPDB, pch = 9, col = "blue3", cex = 2)
      points(c(6, 6), Camels$d18O.VPDB, pch = 4, col = "blue3", cex = 2)
      points(rep(7, length = 9), Horses$d18O.VPDB, pch = 17, col = "blue3", cex = 2)
      legend("topleft", legend = c(expression(bolditalic("Teleoceras major")), expression(bold("Serially sampled")), "UNSM 52272 'Woofy'", "UNSM 27805 'Grandpa John'", "UNSM 52286 'Mary Anning'", "UNSM 27807 'Sparky'", NA, expression(bold("Bulk sampled")),"Males (N = 5)", "Females (N = 8)", NA, NA), col = c(NA, NA, "white", "blue3", "white", "blue3", NA, NA, "black", "black", NA, NA), pt.bg = c(NA, NA, "blue3", "white", "blue3", "white", NA, NA, "blue3", "blue3", NA, NA), pch = c(NA, NA, 22, 22, 21, 21, NA, NA, 22, 21, NA, NA), cex = 1, ncol = 2, pt.cex = 2, text.width = c(2.25, 1.25))
      legend("bottomright", legend = c(expression(bold("Co-occurring ungulates")), "Blastomerycinae (N = 3)", "Camelidae (N = 2)", "Equidae (N = 9)"), col = c(NA, "blue3", "blue3", "blue3"), pch = c(NA, 9, 4, 17), cex = 1, ncol = 1, pt.cex = 2)
    }
    # 87Sr/86Sr
    {
      boxplot(m2.serial$Sr.iso, m3.serial$Sr.iso, m2.bulk$Sr.iso, m3.bulk$Sr.iso, Longi$Sr.iso, Camels$Sr.iso, Horses$Sr.iso, range = 0, ylim = c(0.70856, 0.70903), xlab = NA, ylab = Sr.lab, col = "white", whisklty = "solid", xaxt = "n", cex.lab = 1.5)
      axis(side = 1, at = c(1:4), labels = as.vector(c(expression("M"[2]*"'s"), expression("M"[3]*"'s"),expression("M"[2]*"'s"),expression("M"[3]*"'s"))))
      arrows(1, min(m2.serial$Sr.iso, na.rm = T) - 0.00003, 1, max(m2.serial$Sr.iso, na.rm = T) + 0.00003, length = 0.05, angle = 90, code = 3, lty = "dotted")
      arrows(2, min(m3.serial$Sr.iso) - 0.00003, 2, max(m3.serial$Sr.iso) + 0.00003, length = 0.05, angle = 90, code = 3, lty = "dotted")
      arrows(3, min(m2.bulk$Sr.iso) - 0.00003, 3, max(m2.bulk$Sr.iso) + 0.00003, length = 0.05, angle = 90, code = 3, lty = "dotted")
      arrows(4, min(m3.bulk$Sr.iso) - 0.00003, 4, max(m3.bulk$Sr.iso) + 0.00003, length = 0.05, angle = 90, code = 3, lty = "dotted")
      arrows(5, min(Longi$Sr.iso) - 0.00003, 5, max(Longi$Sr.iso) + 0.00003, length = 0.05, angle = 90, code = 3, lty = "dotted")
      arrows(6, min(Camels$Sr.iso) - 0.00003, 6, max(Camels$Sr.iso) + 0.00003, length = 0.05, angle = 90, code = 3, lty = "dotted")
      arrows(7, min(Horses$Sr.iso) - 0.00003, 7, max(Horses$Sr.iso) + 0.00003, length = 0.05, angle = 90, code = 3, lty = "dotted")
      points(Woofy.results$x, Woofy.results$Sr.iso, pch = 22, bg = "darkorange", col = "white", cex = 2)
      points(John.results$x, John.results$Sr.iso, pch = 22, bg = "white", col = "darkorange", cex = 2)
      points(Mary.results$x, Mary.results$Sr.iso, pch = 21, bg = "darkorange", col = "white", cex = 2)
      points(Sparky.results$x, Sparky.results$Sr.iso, pch = 21, bg = "white", col = "darkorange", cex = 2)
      points(male.bulk$x, male.bulk$Sr.iso, pch = 22, bg = "darkorange", col = "black", cex = 2)
      points(female.bulk$x, female.bulk$Sr.iso, pch = 21, bg = "darkorange", col = "black", cex = 2)
      points(c(5, 5, 5), Longi$Sr.iso, pch = 9, col = "darkorange", cex = 2)
      points(c(6, 6), Camels$Sr.iso, pch = 4, col = "darkorange", cex = 2)
      points(rep(7, length = 9), Horses$Sr.iso, pch = 17, col = "darkorange", cex = 2)
      legend("topleft", legend = c(expression(bolditalic("Teleoceras major")), expression(bold("Serially sampled")), "UNSM 52272 'Woofy'", "UNSM 27805 'Grandpa John'", "UNSM 52286 'Mary Anning'", "UNSM 27807 'Sparky'", NA, expression(bold("Bulk sampled")),"Males (N = 5)", "Females (N = 8)", NA, NA), col = c(NA, NA, "white", "darkorange", "white", "darkorange", NA, NA, "black", "black", NA, NA), pt.bg = c(NA, NA, "darkorange", "white", "darkorange", "white", NA, NA, "darkorange", "darkorange", NA, NA), pch = c(NA, NA, 22, 22, 21, 21, NA, NA, 22, 21, NA, NA), cex = 1, ncol = 2, pt.cex = 2, text.width = c(2.25, 1.25))
      legend("bottomright", legend = c(expression(bold("Co-occurring ungulates")), "Blastomerycinae (N = 3)", "Camelidae (N = 2)", "Equidae (N = 9)"), col = c(NA, "darkorange", "darkorange", "darkorange"), pch = c(NA, 9, 4, 17), cex = 1, ncol = 1, pt.cex = 2)
    }
    dev.off()
  }
