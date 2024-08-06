IE_data <- read_xlsx("all-ie.xlsx")

###############################################################################################################
############################################ Plot for Malaysia ################################################
###############################################################################################################
ie.mal <- IE_data[IE_data$Country == "Malaysia",]

# Getting unique periods and creating an index
ordered_periods <- factor(unique(ie.mal$Period), levels = c("1-6", "7-12", "13-18", "19-24"))

# Original paper estimates 
pooled_estimate <- 40.0
u_pooled_estimate <- 64.59
l_pooled_estimate <- 5.06
pooled_estimate_position <- max(period_index) + 1

# Offset for upper limit of the x-axis
num_breaks <- length(unique_periods) + 1
breaks <- seq(1, num_breaks)
x_range <- c(1, length(levels(ordered_periods)) + 1)
y_range <- c(0, 100)

# Creating an empty plot with adjusted xlim and ylim
par(mar = c(5.1, 5.1, 3.1, 7.1))

plot(x = NULL, y = NULL, type = "n", 
     xlab = "Relative Period After Intervention (Months)", ylab = "Intervention Efficacy (%)", 
     cex.lab = 1.8,
     xlim = x_range, ylim = y_range, 
     xaxt = "n", yaxt = "n")

# Customise x-axis with period labels 
par(cex.axis = 1.7)
axis(1, at = 1:length(levels(ordered_periods)) + 0.2, labels = levels(ordered_periods), las = 1)
colors <- setNames(c("red", "orange", "darkgreen", "darkblue", "purple" ), unique(ie.mal$Method))

for (i in seq_along(unique(ie.mal$Method))) {
  subset_data <- ie.mal[ie.mal$Method == unique(ie.mal$Method)[i], ]
  fixed_period <- match(subset_data$Period, unique_periods) + (i - 1) * 0.1
  
  for (j in seq_along(fixed_period)) {
    points(fixed_period[j], subset_data$IE[j], col = colors[i], pch = 19, cex = 1.1)
    u.midpoint <- ifelse(subset_data$u.IE[j] >= 100, 100, subset_data$u.IE[j])
    l.midpoint <- ifelse(subset_data$l.IE[j] <= 20, 20, subset_data$l.IE[j])
    
    midpoint <- (u.midpoint  + l.midpoint)/2
    
    # Lower bound 
    arrows(x0=fixed_period[j], 
           y0=midpoint, 
           x1=fixed_period[j], 
           y1=ifelse(subset_data$l.IE[j] <= 0, 0, subset_data$l.IE[j]), 
           code=2, 
           angle=ifelse(subset_data$l.IE[j] <= 0, 45, 90), 
           length = 0.05, 
           col=colors[i], 
           lwd=1)
    
    # Upper bound
    arrows(x0=fixed_period[j], 
           y0= midpoint, 
           x1=fixed_period[j], 
           y1=ifelse(subset_data$u.IE[j] >= 100, 102.5, subset_data$u.IE[j]), 
           code=2, 
           angle= ifelse(subset_data$u.IE[j] >= 100, 45, 90),
           length = 0.05, 
           col=colors[i], 
           lwd=1)
    
    if (subset_data$u.IE[j] > 100) {
      text(x = fixed_period[j], 
           y = 102, 
           labels = formatC(subset_data$u.IE[j], format = "f", digits = 2), 
           cex = 1.0,  
           pos = 4)  
    }
  }
}

points(5, pooled_estimate, col = "black", pch = 19, cex = 1.1)

arrows(x0=5,
       y0=l_pooled_estimate,
       x1=5,
       y1=u_pooled_estimate,
       code=3,
       angle=90,
       length = 0.05,
       col="black",
       lwd=1)

labels <- c(levels(unique_periods), "Original Estimate")
axis(1, at = 5, labels = "Original", las = 1)

# Adding legend
#legend(x = "right", y = 0.5, legend = names(colors), fill = colors, title = "Method", xpd = TRUE, inset=c(-0.15,0), cex = 0.95)

abline(v = seq(1.73, max(breaks), by = 1), lty = 2)

par(las=1)
y_labels <- seq(0, 100, by = 10)
axis(2, at = y_labels, labels = y_labels, cex = 1.7)

title(main = "Malaysia", cex.main = 2.0)

###############################################################################################################
############################################### Plot for Brazil ###############################################
###############################################################################################################
ie.brazil <- IE_data[IE_data$Country == "Brazil",]

# Getting unique periods and creating an index
ordered_periods <- factor(unique(ie.brazil$Period), levels = c("1-6", "7-12", "13-18", "19-24"))

# Creating an empty plot 
label_offset <- (length(unique(ie.brazil$Method)) - 1) * 0.1 / 2
additional_offset <- label_offset

# Original paper estimates 
pooled_estimate <- 69.4
u_pooled_estimate <- 79.40
l_pooled_estimate <- 54.4
pooled_estimate_position <- max(period_index) + 1

# Offset for upper limit of the x-axis
num_breaks <- length(unique_periods) + 1
breaks <- seq(1, num_breaks)
x_range <- c(1, length(levels(ordered_periods)) + 1)
y_range <- c(-50, 100)

# Creating an empty plot with adjusted xlim and ylim
par(mar = c(5.1, 5.1, 3.1, 7.1))

plot(x = NULL, y = NULL, type = "n", 
     xlab = "Relative Period After Intervention (Months)", ylab = "Intervention Efficacy (%)", 
     cex.lab = 1.8,
     xlim = x_range, ylim = y_range, 
     xaxt = "n", yaxt = "n")

# Customise x-axis with period labels 
par(cex.axis = 1.7)
axis(1, at = 1:length(levels(ordered_periods)) + 0.2, labels = levels(ordered_periods), las = 1)

colors <- setNames(c("red", "orange", "darkgreen", "darkblue", "purple" ), unique(ie.brazil$Method))

for (i in seq_along(unique(ie.brazil$Method))) {
  subset_data <- ie.brazil[ie.brazil$Method == unique(ie.brazil$Method)[i], ]
  fixed_period <- match(subset_data$Period, unique_periods) + (i - 1) * 0.1
  
  for (j in seq_along(fixed_period)) {
  
  if (subset_data$IE[j] < -50) {
    arrows(x0=fixed_period[j], 
           y0=-49, 
           x1=fixed_period[j], 
           y1=-50, 
           col=colors[i], 
           angle = 45,
           length = 0.05,  # adjust this value to make the arrow smaller or larger
           lwd=4)
    text(x = fixed_period[j], 
         y = -50,  # adjust as needed
         labels = formatC(subset_data$IE[j], format = "f", digits = 2), 
         cex = 1.4, 
         pos = ifelse(i == 2, 4, 4))  # places the text below
  } else {
    points(fixed_period[j], subset_data$IE[j], col = colors[i], pch = 19, cex = 1.1)
  }
  
  u.midpoint <- ifelse(subset_data$u.IE[j] >= 100, 100, subset_data$u.IE[j])
  l.midpoint <- ifelse(subset_data$l.IE[j] <= -50, -50, subset_data$l.IE[j])
  
  midpoint <- (u.midpoint + l.midpoint) / 2
  
  # Lower bound 
  arrows(x0=fixed_period[j], 
         y0=midpoint, 
         x1=fixed_period[j], 
         y1=ifelse(subset_data$l.IE[j] <= -50, -55, subset_data$l.IE[j]), 
         code=2, 
         angle=ifelse(subset_data$l.IE[j] <= -50, -90, 90), 
         length = 0.05, 
         col=colors[i], 
         lwd=1)
  
  # Upper bound
  arrows(x0=fixed_period[j], 
         y0= midpoint, 
         x1=fixed_period[j], 
         y1=ifelse(subset_data$u.IE[j] >= 100, 102.5, subset_data$u.IE[j]), 
         code=2, 
         angle= ifelse(subset_data$u.IE[j] >= 100, 45, 90),
         length = 0.05, 
         col=colors[i], 
         lwd=1)
  
  if (subset_data$u.IE[j] > 100) {
    text(x = fixed_period[j], 
         y = 100,
         labels = formatC(subset_data$u.IE[j], format = "f", digits = 2), 
         cex = 1.4,  
         pos = 4)  
  }
  
  if (subset_data$l.IE[j] < -50) {
    text(x = fixed_period[j], 
         y = -55,  
         labels = formatC(subset_data$l.IE[j], format = "f", digits = 2), 
         cex = 1.4,  
         pos = 4) 
  }
  }
}

points(5, pooled_estimate, col = "black", pch = 19, cex = 1.1)

arrows(x0=5,
       y0=l_pooled_estimate,
       x1=5,
       y1=u_pooled_estimate,
       code=3,
       angle=90,
       length = 0.05,
       col="black",
       lwd=1)

labels <- c(levels(unique_periods), "Original Estimate")
axis(1, at = 5, labels = "Original", las = 1)

# Adding legend
#legend(x = "right", y = 0.5, legend = names(colors), fill = colors, title = "Method", xpd = TRUE, inset=c(-0.18,0), cex = 0.7)

abline(h = 0, v = seq(1.73, max(breaks), by = 1), lty = 2, lwd = 1)

par(las=1)
y_labels <- seq(-50, 100, by = 10)
axis(2, at = y_labels, labels = y_labels)

title(main = "Brazil", cex.main = 2.0)

###############################################################################################################
############################################# Plot for Singapore ##############################################
###############################################################################################################
ie.sg <- IE_data[IE_data$Country == "Singapore",]

# Getting unique periods and creating an index
ordered_periods <- factor(unique(ie.sg$Period), levels = c("1-6", "7-12", "13-18", "19-24"))
# Creating an empty plot 
label_offset <- (length(unique(ie.mal$Method)) - 1) * 0.1 / 2
additional_offset <- label_offset

# Original paper estimates 
pooled_estimate <- 56.88
u_pooled_estimate <- 58.46
l_pooled_estimate <- 51.88
pooled_estimate_position <- max(period_index) + 1

# Offset for upper limit of the x-axis
num_breaks <- length(unique_periods) + 1
breaks <- seq(1, num_breaks)
x_range <- c(1, length(levels(ordered_periods)) + 1)
y_range <- c(-80, 100)

# Creating an empty plot with adjusted xlim and ylim
plot(x = NULL, y = NULL, type = "n", 
     xlab = "Relative Period After Intervention (Months)", ylab = "Intervention Efficacy (%)", 
     cex.lab = 1.8,
     xlim = x_range, ylim = y_range, 
     xaxt = "n", yaxt = "n")

# Customise x-axis with period labels 
par(cex.axis = 1.7)
axis(1, at = 1:length(levels(ordered_periods)) + 0.2, labels = levels(ordered_periods), las = 1)

colors <- setNames(c("red", "orange", "darkgreen", "darkblue", "purple" ), unique(ie.mal$Method))

for (i in seq_along(unique(ie.sg$Method))) {
  subset_data <- ie.sg[ie.sg$Method == unique(ie.sg$Method)[i], ]
  fixed_period <- match(subset_data$Period, unique_periods) + (i - 1) * 0.1
  
  for (j in seq_along(fixed_period)) {
    if (subset_data$IE[j] < -80) {
      arrows(x0=fixed_period[j], 
             y0=-72, 
             x1=fixed_period[j], 
             y1=-73, 
             col=colors[i], 
             angle = 45,
             length = 0.05, 
             lwd=4)
      text(x = fixed_period[j], 
           y = -73, 
           labels = formatC(subset_data$IE[j], format = "f", digits = 2), 
           cex = 1.4, 
           pos = 4)  
    } else {
      points(fixed_period[j], subset_data$IE[j], col = colors[i], pch = 19, cex = 1.1)
    }
    
    u.midpoint <- ifelse(subset_data$u.IE[j] >= 100, 100, subset_data$u.IE[j])
    l.midpoint <- ifelse(subset_data$l.IE[j] <= -80, -80, subset_data$l.IE[j])
    
    midpoint <- (u.midpoint + l.midpoint) / 2
    
    # Lower bound 
    arrows(x0=fixed_period[j], 
           y0=midpoint, 
           x1=fixed_period[j], 
           y1=ifelse(subset_data$l.IE[j] <= -80, -85, subset_data$l.IE[j]), 
           code=2, 
           angle=ifelse(subset_data$l.IE[j] <= -80, 45, 90), 
           length = 0.05, 
           col=colors[i], 
           lwd=1)
    
    if (subset_data$l.IE[j] < -80) {
      text(x = fixed_period[j], 
           y = -83, 
           labels = formatC(subset_data$l.IE[j], format = "f", digits = 2), 
           cex = 1.4, 
           pos = 4)  
    }
    
    # Upper bound
    arrows(x0=fixed_period[j], 
           y0= midpoint, 
           x1=fixed_period[j], 
           y1=ifelse(subset_data$u.IE[j] >= 100, 102.5, subset_data$u.IE[j]), 
           code=2, 
           angle= ifelse(subset_data$u.IE[j] >= 100, 45, 90),
           length = 0.05, 
           col=colors[i], 
           lwd=1)
    
    if (subset_data$u.IE[j] > 100) {
      text(x = fixed_period[j], 
           y = 102.5, 
           labels = formatC(subset_data$u.IE[j], format = "f", digits = 2), 
           cex = 1.4, 
           pos = 4)  
    }
  }
}

points(5, pooled_estimate, col = "black", pch = 19, cex = 1.1)

arrows(x0=5,
       y0=l_pooled_estimate,
       x1=5,
       y1=u_pooled_estimate,
       code=3,
       angle=90,
       length = 0.05,
       col="black",
       lwd=1)

labels <- c(levels(unique_periods), "Original Estimate")
axis(1, at = 5, labels = "Original", las = 1)

# Adding legend
#legend(x = "right", y = 0.5, legend = names(colors), fill = colors, title = "Method", xpd = TRUE, inset=c(-0.18,0), cex = 0.8)

abline(h = 0, v = seq(1.73, max(breaks), by = 1), lty = 2, lwd = 1)

par(las=1)
y_labels <- seq(-80, 100, by = 10)
axis(2, at = y_labels, labels = y_labels, cex = 1.7)

title(main = "Singapore", cex.main = 2.0)



