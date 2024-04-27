library(XICOR)
library(FNN)
library(MASS)
library(ggplot2)
library(gridExtra)

mic <- function(x,y,k = 2){
  return(sqrt(1-exp(-2*max(0,mutinfo(x,y,k)))))
}


# Function to generate donut data
gen_donut <- function(n, th) {
  outer_radius <- 1   # Outer radius of the donut
  inner_radius <- 1 - th  # Inner radius of the donut
  
  # Generate random angles
  angles <- runif(n, min = 0, max = 2 * pi)
  
  # Generate random radii
  radii <- sqrt(runif(n, min = inner_radius^2, max = outer_radius^2))
  
  # Convert polar coordinates to Cartesian coordinates
  x <- radii * cos(angles)
  y <- radii * sin(angles)
  return(data.frame(x = x, y = y))
}

# Generate data
set.seed(1000)
d1 <- gen_donut(10000, 1)
d2 <- gen_donut(10000, 0.5)
d3 <- gen_donut(10000, 0)

# Calculate statistics
stats1 <- list(cor = cor(d1$x, d1$y), mic = mic(d1$x, d1$y), xicor = xicor(d1$x, d1$y))
stats2 <- list(cor = cor(d2$x, d2$y), mic = mic(d2$x, d2$y), xicor = xicor(d2$x, d2$y))
stats3 <- list(cor = cor(d3$x, d3$y), mic = mic(d3$x, d3$y), xicor = xicor(d3$x, d3$y))

# Define a plotting function
plot_donut <- function(data, stats) {
  p <- ggplot(data, aes(x = x, y = y)) +
    geom_point(alpha = 0.4, color = "purple") +
    ggtitle(sprintf("??: %.4f,\n R: %.4f,\n ??: %.4f", stats$cor, stats$mic, stats$xicor)) +
    theme_minimal() +
    theme(plot.title = element_text(size = 28, hjust = 0.5),
          axis.title.x = element_blank(),  # Remove x-axis title
          axis.title.y = element_blank(),  # Remove y-axis title
          axis.text.x = element_blank(),   # Remove x-axis text (labels)
          axis.text.y = element_blank(),   # Remove y-axis text (labels)
          axis.ticks = element_blank())
  return(p)
}

# Create plots
p1 <- plot_donut(d1, stats1)
p2 <- plot_donut(d2, stats2)
p3 <- plot_donut(d3, stats3)

# Arrange plots in a row
grid_plot <- grid.arrange(p1, p2, p3, nrow = 1)

# Save the arranged plot to a PNG file
ggsave("combined_plot.png", plot = grid_plot, width = 16, height = 6.5, dpi = 300)

################################################################################
  
gen_parab <- function(n,noisel){
  x <- runif(n, -1, 1)
  y <- x^2 + runif(n, -noisel, noisel)
  return(data.frame(x = x, y = y))
}

# Generate data
set.seed(2000)
d1 <- gen_parab(10000, 1.5)
d2 <- gen_parab(10000, 0.75)
d3 <- gen_parab(10000, 0)

# Calculate statistics
stats1 <- list(cor = cor(d1$x, d1$y), mic = mic(d1$x, d1$y), xicor = xicor(d1$x, d1$y),
               cori = cor(d1$y, d1$x), mici = mic(d1$y, d1$x), xicori = xicor(d1$y, d1$x))
stats2 <- list(cor = cor(d2$x, d2$y), mic = mic(d2$x, d2$y), xicor = xicor(d2$x, d2$y),
               cori = cor(d2$y, d2$x), mici = mic(d2$y, d2$x), xicori = xicor(d2$y, d2$x))
stats3 <- list(cor = cor(d3$x, d3$y), mic = mic(d3$x, d3$y), xicor = xicor(d3$x, d3$y),
               cori = cor(d3$y, d3$x), mici = mic(d3$y, d3$x), xicori = xicor(d3$y, d3$x))


# Define a plotting function
plot_parab <- function(data, stats) {
  p <- ggplot(data, aes(x = x, y = y)) +
    geom_point(alpha = 0.4, color = "orange") +
    ggtitle(sprintf("??: %.4f,\n R: %.4f,\n ??(x,y): %.4f,\n ??(y,x): %.4f", stats$cor, stats$mic, stats$xicor, stats$xicori)) +
    theme_minimal() +
    theme(plot.title = element_text(size = 28, hjust = 0.5),
          axis.title.x = element_blank(),  # Remove x-axis title
          axis.title.y = element_blank(),  # Remove y-axis title
          axis.text.x = element_blank(),   # Remove x-axis text (labels)
          axis.text.y = element_blank(),   # Remove y-axis text (labels)
          axis.ticks = element_blank())
  return(p)
}

# Create plots
p1 <- plot_parab(d1, stats1)
p2 <- plot_parab(d2, stats2)
p3 <- plot_parab(d3, stats3)

# Arrange plots in a row
grid_plot <- grid.arrange(p1, p2, p3, nrow = 1)

# Save the arranged plot to a PNG file
ggsave("combined_plot.png", plot = grid_plot, width = 16, height = 6.5, dpi = 300)

################################################################################

gen_norm <- function(n, rho){
  # Define the number of samples, means, and desired correlation
  means <- c(0, 0)  # Means of the two variables
  
  # Create the covariance matrix based on the specified correlation
  cov_matrix <- matrix(c(1, rho, rho, 1), nrow=2)
  
  # Generate the correlated normal variables
  correlated_data <- as.data.frame(mvrnorm(n=n, mu=means, Sigma=cov_matrix))
  colnames(correlated_data) <- c("x", "y")
  return(correlated_data)
}

# Generate data
set.seed(3000)
d1 <- gen_norm(10000, 0.4)
d2 <- gen_norm(10000, 0.7)
d3 <- gen_norm(10000, 1)

# Calculate statistics
stats1 <- list(cor = cor(d1$x, d1$y), mic = mic(d1$x, d1$y), xicor = xicor(d1$x, d1$y))
stats2 <- list(cor = cor(d2$x, d2$y), mic = mic(d2$x, d2$y), xicor = xicor(d2$x, d2$y))
stats3 <- list(cor = cor(d3$x, d3$y), mic = mic(d3$x, d3$y), xicor = xicor(d3$x, d3$y))

# Define a plotting function
plot_norm <- function(data, stats) {
  p <- ggplot(data, aes(x = x, y = y)) +
    geom_point(alpha = 0.8, color = "lightblue") +
    ggtitle(sprintf("??: %.4f,\n R: %.4f,\n ??: %.4f", stats$cor, stats$mic, stats$xicor)) +
    theme_minimal() +
    theme(plot.title = element_text(size = 28, hjust = 0.5),
          axis.title.x = element_blank(),  # Remove x-axis title
          axis.title.y = element_blank(),  # Remove y-axis title
          axis.text.x = element_blank(),   # Remove x-axis text (labels)
          axis.text.y = element_blank(),   # Remove y-axis text (labels)
          axis.ticks = element_blank())
  return(p)
}

# Create plots
p1 <- plot_norm(d1, stats1)
p2 <- plot_norm(d2, stats2)
p3 <- plot_norm(d3, stats3)

# Arrange plots in a row
grid_plot <- grid.arrange(p1, p2, p3, nrow = 1)

# Save the arranged plot to a PNG file
ggsave("combined_plot.png", plot = grid_plot, width = 16, height = 6.5, dpi = 300)


################################################################################

library(microbenchmark)

gen_unif <- function(n){
  x <- runif(n, 0, 1)
  y <- runif(n, 0, 1)
  return(data.frame(x = x, y = y))
}


n_values <- c(100, 1000, 10000, 50000)

# Initialize a list to store results
resc <- resm <- resx <- matrix(ncol = 4, nrow = 6)

# Loop over the range of n values
for (i in 1:length(n_values)) {
  print(i)
  n = n_values[i]

  for(j in 1:6){
    d <- gen_unif(n)
    resc[j,i] <- microbenchmark(cor(d$x, d$y), times = 1)[1,2]/10^6
    resm[j,i] <- microbenchmark(mic(d$x, d$y), times = 1)[1,2]/10^6
    resx[j,i] <- microbenchmark(xicor(d$x, d$y), times = 1)[1,2]/10^6
  }

}


# Create data frames

resc <- data.frame(x = n_values, y = apply(resc[2:6,], 2, mean))
resm <- data.frame(x = n_values, y = apply(resm[2:6,], 2, mean))
resx <- data.frame(x = n_values, y = apply(resx[2:6,], 2, mean))

# Plotting function using ggplot2
plot_data <- function(data, title) {
  ggplot(data, aes(x = x, y = y)) +
    geom_point(size=4, shape=19, colour="purple") +  # Larger, blue points
    geom_line(size=1.5, colour="lightblue") +  # Thicker, red line
    geom_point() +
    geom_line() +
    scale_x_log10() +
    scale_y_log10() +
    labs(title = title, x = "Points", y = "Time, ms") +
    theme(plot.title = element_text(size = 28, hjust = 0.5),
          axis.title.x = element_text(size = 24, vjust = 0.5), 
          axis.title.y = element_text(size = 24, vjust = 0.5), 
          axis.text.x = element_text(size = 24),
          axis.text.y = element_text(size = 24))
}


# Create plots
p1 <- plot_data(resc, "??")
p2 <- plot_data(resm, "R")
p3 <- plot_data(resx, "??")

# Arrange plots in a row
grid_plot <- grid.arrange(p1, p2, p3, nrow = 1)

ggsave("combined_plot.png", plot = grid_plot, width = 16, height = 4.5, dpi = 300)


################################################################################


resc <- resm <- resx <- matrix(ncol = 4, nrow = 3)
npts <- c(100,1000,10000)
rhos <- c(0.3, 0.5, 0.7, 0.9)

for(i in 1:3){
  for(j in 1:4){
    print(j)
    cors <- mics <- xicors <- c()
    for(k in 1:25){
      d <- gen_norm(npts[i], rhos[i])
      cors <- c(cors, cor(d[,1],d[,2]))
      mics <- c(mics, mic(d[,1],d[,2], k = 2))
      xicors <- c(xicors, xicor(d[,1],d[,2]))
    }
    resc[i,j] <- sd(cors)/mean(cors)
    resm[i,j] <- sd(mics)/mean(mics)
    resx[i,j] <- sd(xicors)/mean(xicors)
  }
}

write.csv(round(resc, 3), "cor.csv")
write.csv(round(resm, 3), "corm.csv")
write.csv(round(resx, 3), "corx.csv")


