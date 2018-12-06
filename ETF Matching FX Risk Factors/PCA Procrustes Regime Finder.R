library(quantmod)
library(sqldf)
library(base)
library(stats)
library(pracma)
library(utils)
library(tcltk)
library(lubridate)

#This script clusters the currency data into a given number of regimes (undetermined length) using subspace angle of the top five
#principal components on rolling window returns as a decision boundary. A "best fit" rolling window is assigned for each regime and
# each of the rolling windows in the regime are rotated towards the best fit using procrustes rotation. This is so that we
# are creating stable time series of the principal component scores to improve the interpretability of the results when
#we correlate them with ETFs.

#Output is the rotated factors for each window, and the rotated factor scores. 

returns2<- read.csv("C:/Users/Shane/Documents/R/FX Data - All crosses.csv",header = TRUE)

x<-parse_date_time(returns2[,1],orders="ymd")
x<-as.Date(x)
#set up as a time series object?
currencies<- xts(returns2[,2:172],order.by=x)

window<-120
length=dim(currencies)[1]-window
factor_scores <- xts(matrix(,nrow=dim(currencies)[1],ncol=5*130), order.by=x)
#write.csv(factor_scores, file="f2s.csv",row.names = TRUE)

comps_rot<-18
threshold <- 0.17

factor_hist <-matrix(,nrow=171,ncol=100*comps_rot)
factors_Rcor <- matrix(,nrow=171,ncol=length)
factors_cor <- matrix(,nrow=171,ncol=length)

best_angle <- matrix(,nrow=10,ncol=length)
current_reg_info <- matrix(,nrow=10,ncol=260)
W1_top5 <- matrix(,nrow=171,ncol=5)
W2_top5 <- matrix(,nrow=171,ncol=5)
W1_rot <- matrix(,nrow=171, ncol=comps_rot)
W1_try <- matrix(,nrow=171, ncol=comps_rot)
W1_cur <- matrix(,nrow=171, ncol=comps_rot)
W2_rot <- matrix(,nrow=171, ncol=comps_rot)

window_info <- matrix(,nrow=23,ncol=200)
disp <- 0
num_regimes <- 0
new_regime <-1
cur_best_window <- 1
CR_start <- 1
CR_length <- 0
CR_pause <- 0
best_try <- 0

Current_index <- 1

currencies2 <-xts()

for (i in 1:130){
	currencies2<-xts()
	for (ii in 1:171)
	{
		currencies2 <- merge(currencies2, currencies[i:(i+window-1),ii])
		ii=ii+1
	}
	corC <- cor(currencies2)
	compC <-svd(corC)
	W2_top5<-compC$u[,1:5]
	W2_Rtop5<-compC$u[,1:5]
	W2_rot <- compC$u[,1:comps_rot]
	A <- compC$u[,1:comps_rot]
	if ( new_regime ==1){
		cur_best_window <- i
		CR_start <- i
		num_regimes <- num_regimes+1
		CR_length <- 1
		W1_top5[1:171,1:5] <- W2_top5[1:171,1:5]
		W1_rot[1:171, 1:comps_rot] <- W2_rot[1:171, 1:comps_rot]
		current_reg_info <- matrix(,nrow=10,ncol=260)
		new_regime <- 0
	} else {
		CR_length =CR_length+1
	}

	At = t(A)
	AtB=At %*% W1_rot
	compC2=svd(AtB)
	W2_rot = A %*% compC2$u %*% t(compC2$v)
	condition = subspace(W1_top5[,1],W2_rot[,1])

	if (condition >= threshold) { #if the angle between current B and i is greater than 0.17 find new B
		angle_reg <- matrix(,nrow=CR_length,ncol=CR_length)
		angle_avg <- matrix(,nrow=CR_length,ncol=2)
		cur_avg <- mean(current_reg_info[6,1:(CR_length-1)])
		best_avg <- cur_avg
		best_try <- 0

		for (j in 0:(CR_length-1))	{

				tw1<-xts()
				for (ii in 1:171)
				{
					tw1 <- merge(tw1, currencies[(CR_start+j):(CR_start+j+window-1),ii])
					ii=ii+1
				}
				cor_tw1 <- corr(tw1)
				tempSVD <- svd(cor_tw1)
				tw1_top5 <- tempSVD$u[1:171,1:5]
				tw1_rot <- tempSVD$u[1:171,1:5]

			for (jj in 0:(CR_length-1))	{
				tw2<-xts()
				for (ii in 1:171)
				{
					tw2 <- merge(tw2, currencies[(CR_start+jj):(CR_start+jj+window-1),ii])
					ii=ii+1
				}
				cor_tw2 <- cor(tw2)
				tempSVD2 <- svd(cor_tw2)
				At <- tempSVD$u[1:171,1:comps_rot]
				tAt <- t(tempSVD$u[1:171,1:comps_rot])
				tAtB <- tAt %*% tw1_rot
				tempSVD3 <- svd(tAtB)
				tw2_rot <- At %*% tempSVD3$u %*% t(tempSVD3$v)
				angle_reg[(j+1),(jj+1)] <- subspace(tw1_top5[1:171,1],tw2_rot[1:171,1])
			}
			angle_avg[(j+1),1] <- mean(angle_reg[(j+1),1:CR_length])
			angle_avg[(j+1),2]<-max(angle_reg[(j+1),1:CR_length])
	}

	for (j in 1:CR_length)	{
		if (angle_avg[j,2] < threshold)	{
			if	(angle_avg[j,2]<threshold)	{
				best_try<-j
				best_avg <- angle_avg[j,1]
			}
		}
	}
	if(best_try>0) {	#adjust best window (B in procrustes rotation) in regime if a new one is found
		cur_best_window <- CR_start+best_try-1
		tw1<-xts()
		for (ii in 1:171)
		{
			tw1 <- merge(tw1, currencies[cur_best_window:(cur_best_window+window-1),ii])
			ii=ii+1
		}
		cor_tw1 <- cor(tw1)
		tempSVD <- svd(cor_tw1)
		W1_top5 <- tempSVD$u[,1:5]
		W1_rot <- tempSVD$u[,1:comps_rot]
		Current_index <- 1
		for (jfk in CR_start:i) {
			tw2<-xts()
			for (ii in 1:171)
			{
				tw2 <- merge(tw2, currencies[jfk:(jfk+window-1),ii])
				ii=ii+1
			}
			cor_tw2 <- cor(tw2)
			tempSVD2 <- svd(cor_tw2)
			tw2_top5 <- tempSVD2$u[,1:5]
			A <- tempSVD2$u[,1:comps_rot]
			At <- t(A)
			AtB <- At %*% W1_rot
			tempSVD3 <- svd(AtB)
			tw2_rot <- A %*% tempSVD3$u %*% t(tempSVD3$v)
			factors_cor <- tw2_top5[,1:5]
			factors_Rcor <- tw2_rot[,1]
			current_reg_info[1,Current_index] <- subspace(W1_top5[,1],tw1_top5[,1])
			current_reg_info[2,Current_index] <- subspace(W1_top5[,2],tw1_top5[,2])
			current_reg_info[3,Current_index] <- subspace(W1_top5[,3],tw1_top5[,3])
			current_reg_info[4,Current_index] <- subspace(W1_top5[,4],tw1_top5[,4])
			current_reg_info[5,Current_index] <- subspace(W1_top5[,5],tw1_top5[,5])
			tw2_top5 <- tw2_rot[,1:5]
			current_reg_info[6,Current_index] <- subspace(W1_top5[,1],tw1_top5[,1])
			current_reg_info[7,Current_index] <- subspace(W1_top5[,2],tw1_top5[,2])
			current_reg_info[8,Current_index] <- subspace(W1_top5[,3],tw1_top5[,3])
			current_reg_info[9,Current_index] <- subspace(W1_top5[,4],tw1_top5[,4])
			current_reg_info[10,Current_index] <- subspace(W1_top5[,5],tw1_top5[,5])
			Current_index <- Current_index +1

		}

	} else { #if no new window was found, close out the regime and save rotated scores for ETF correlation
		best_angle[1:10,CR_start:(i-1)] <- current_reg_info[1:10,1:(CR_length-1)]
		window_info[1,num_regimes] <- mean(current_reg_info[1,1:(CR_length-1)])
		window_info[2,num_regimes] <- mean(current_reg_info[2,1:(CR_length-1)])
		window_info[3,num_regimes] <- mean(current_reg_info[3,1:(CR_length-1)])
		window_info[4,num_regimes] <- mean(current_reg_info[4,1:(CR_length-1)])
		window_info[5,num_regimes] <- mean(current_reg_info[5,1:(CR_length-1)])
		window_info[6,num_regimes] <- mean(current_reg_info[6,1:(CR_length-1)])
		window_info[7,num_regimes] <- mean(current_reg_info[7,1:(CR_length-1)])
		window_info[8,num_regimes] <- mean(current_reg_info[8,1:(CR_length-1)])
		window_info[9,num_regimes] <- mean(current_reg_info[9,1:(CR_length-1)])
		window_info[10,num_regimes] <- mean(current_reg_info[10,1:(CR_length-1)])
		window_info[11,num_regimes] <- sd(current_reg_info[1,1:(CR_length-1)], na.rm = TRUE)
		window_info[12,num_regimes] <- sd(current_reg_info[2,1:(CR_length-1)], na.rm = TRUE)
		window_info[13,num_regimes] <- sd(current_reg_info[3,1:(CR_length-1)], na.rm = TRUE)
		window_info[14,num_regimes] <- sd(current_reg_info[4,1:(CR_length-1)], na.rm = TRUE)
		window_info[15,num_regimes] <- sd(current_reg_info[5,1:(CR_length-1)], na.rm = TRUE)
		window_info[16,num_regimes] <- sd(current_reg_info[6,1:(CR_length-1)], na.rm = TRUE)
		window_info[17,num_regimes] <- sd(current_reg_info[7,1:(CR_length-1)], na.rm = TRUE)
		window_info[18,num_regimes] <- sd(current_reg_info[8,1:(CR_length-1)], na.rm = TRUE)
		window_info[19,num_regimes] <- sd(current_reg_info[9,1:(CR_length-1)], na.rm = TRUE)
		window_info[20,num_regimes] <- sd(current_reg_info[10,1:(CR_length-1)], na.rm = TRUE)
		window_info[21,num_regimes] <- CR_start
		window_info[22,num_regimes] <- i-1
		window_info[23,num_regimes] <- cur_best_window
		for (jfk in CR_start:(i-1)) {
			tw2<-xts()
			for (ii in 1:171)
			{
				tw2 <- merge(tw2, currencies[jfk:(jfk+window-1),ii])
				ii=ii+1
			}
			cor_tw2 <- cor(tw2)
			tempSVD2 <- svd(cor_tw2)
			tw2_top5 <- tempSVD2$u[,1:5]
			A <- tempSVD2$u[,1:comps_rot]
			At <- t(A)
			AtB <- At %*% W1_rot
			tempSVD3 <- svd(AtB)
			tw2_rot <- A %*% tempSVD3$u %*% t(tempSVD3$v)
			#factors_cor <= tw2_top5[,1:5] #dont want to waste time with unrotated scores
			tw2_top5 <- tw2_rot[,1:5]
			scorez <- tw2 %*% tw2_top5
			for (kk in 1:5)
			{
				factor_scores[jfk: (jfk +window-1),((jfk-1)*5+kk)] <- scorez[1:120,kk]
			}
		}
		cur_best_window <- i
		CR_start <- i
		num_regimes <- num_regimes +1
		CR_length <- 1
		W1_top5 <- compC$u[,1:5] # here we make regimes rather distinct because we start with unrotated factors
		#we could try to add some continuity to regimes if we rotated towards the best window in the prior
		#regime, but this would likely increase the number of regimes
		W1_rot <- compC$u[,1:comps_rot]
		factors_cor[,i] <- W1_top5[,1]
		factors_Rcor[,i] <- W1_rot[,1]
		current_reg_info <- matrix(,nrow=10,ncol=260)
		new_regime <- 0
	}
	} else { #if threshold is not exceed this just takes some stats about the rotation
		current_reg_info[1,CR_length] <- subspace(W1_top5[,1],W2_top5[,1])
		current_reg_info[2,CR_length] <- subspace(W1_top5[,2],W2_top5[,2])
		current_reg_info[3,CR_length] <- subspace(W1_top5[,3],W2_top5[,3])
		current_reg_info[4,CR_length] <- subspace(W1_top5[,4],W2_top5[,4])
		current_reg_info[5,CR_length] <- subspace(W1_top5[,5],W2_top5[,5])
		current_reg_info[6,CR_length] <- subspace(W1_top5[,1],W2_rot[,1])
		current_reg_info[7,CR_length] <- subspace(W1_top5[,2],W2_rot[,2])
		current_reg_info[8,CR_length] <- subspace(W1_top5[,3],W2_rot[,3])
		current_reg_info[9,CR_length] <- subspace(W1_top5[,4],W2_rot[,4])
		current_reg_info[10,CR_length] <- subspace(W1_top5[,5],W2_rot[,5])
	}
}

#need to find the best window in the current regime if threshold wasn't flagged for i = length

angle_reg <- matrix(,nrow=CR_length,ncol=CR_length)
angle_avg <- matrix(,nrow=CR_length,ncol=2)
cur_avg <- mean(current_reg_info[6,1:(CR_length-1)])
best_avg <- cur_avg
best_try <- 0

for (j in 0:(CR_length-1))	{
	tw1<-xts()
	for (ii in 1:171)
	{
		tw1 <- merge(tw1, currencies[(CR_start+j):(CR_start+j+window-1),ii])
		ii=ii+1
	}
	cor_tw1 <- cor(tw1)
	tempSVD <- svd(cor_tw1)
	tw1_top5 <- tempSVD$u[1:171,1:5]
	tw1_rot <- tempSVD$u[1:171,1:5]

	for (jj in 0:(CR_length-1))	{
		tw2<-xts()
		for (ii in 1:171)
		{
			tw2 <- merge(tw2, currencies[(CR_start+jj):(CR_start+jj+window-1),ii])
			ii=ii+1
		}
		cor_tw2 <- cor(tw2)
		tempSVD2 <- svd(cor_tw2)
		At <- tempSVD$u[1:171,1:comps_rot]
		tAt <- t(tempSVD$u[1:171,1:comps_rot])
		tAtB <- tAt %*% tw1_rot
		tempSVD3 <- svd(tAtB)
		tw2_rot <- At %*% tempSVD3$u %*% t(tempSVD3$v)
		angle_reg[(j+1),(jj+1)] <- subspace(tw1_top5[1:171,1],tw2_rot[1:171,1])
	}
	angle_avg[(j+1),1] <- mean(angle_reg[(j+1),1:CR_length])
	angle_avg[(j+1),2]<-max(angle_reg[(j+1),1:CR_length])
	}

for (j in 1:CR_length)
{
		if (angle_avg[j,2] < threshold)	{
			if	(angle_avg[j,2]<threshold)	{
				best_try<-j;
				best_avg <- angle_avg[j,1]
			}
		}
}
if(best_try>0)
{	#adjust best window (B in procrustes rotation) in regime if a new one is found
	cur_best_window <- CR_start+best_try-1
	tw1<-xts()
	for (ii in 1:171)
	{
		tw1 <- merge(tw1, currencies[cur_best_window:(cur_best_window+window-1),ii])
		ii=ii+1
	}
	cor_tw1 <- cor(tw1)
	tempSVD <- svd(cor_tw1)
	W1_top5 <- tempSVD$u[,1:5]
	W1_rot <- tempSVD$u[,1:comps_rot]
	Current_index <- 1
	for (jfk in CR_start:i)
	{
			tw2<-xts()
			for (ii in 1:171)
			{
				tw2 <- merge(tw2, currencies[jfk:(jfk+window-1),ii])
				ii=ii+1
			}
			cor_tw2 <- cor(tw2)
			tempSVD2 <- svd(cor_tw2)
			tw2_top5 <- tempSVD2$u[,1:5]
			A <- tempSVD2$u[,1:comps_rot]
			At <- t(A)
			AtB <- At %*% W1_rot
			tempSVD3 <- svd(AtB)
			tw2_rot <- A %*% tempSVD3$u %*% t(tempSVD3$v)
			factors_cor <- tw2_top5[,1:5]
			factors_Rcor <- tw2_rot[,1]
			current_reg_info[1,Current_index] <- subspace(W1_top5[,1],tw1_top5[,1])
			current_reg_info[2,Current_index] <- subspace(W1_top5[,2],tw1_top5[,2])
			current_reg_info[3,Current_index] <- subspace(W1_top5[,3],tw1_top5[,3])
			current_reg_info[4,Current_index] <- subspace(W1_top5[,4],tw1_top5[,4])
			current_reg_info[5,Current_index] <- subspace(W1_top5[,5],tw1_top5[,5])
			tw2_top5 <- tw2_rot[,1:5]
			current_reg_info[6,Current_index] <- subspace(W1_top5[,1],tw1_top5[,1])
			current_reg_info[7,Current_index] <- subspace(W1_top5[,2],tw1_top5[,2])
			current_reg_info[8,Current_index] <- subspace(W1_top5[,3],tw1_top5[,3])
			current_reg_info[9,Current_index] <- subspace(W1_top5[,4],tw1_top5[,4])
			current_reg_info[10,Current_index] <- subspace(W1_top5[,5],tw1_top5[,5])
			Current_index <- Current_index +1

	}

} #close out the regime and save rotated scores for ETF correlation
best_angle[1:10,CR_start:(i-1)] <- current_reg_info[1:10,1:(CR_length-1)]
window_info[1,num_regimes] <- mean(current_reg_info[1,1:(CR_length-1)])
window_info[2,num_regimes] <- mean(current_reg_info[2,1:(CR_length-1)])
window_info[3,num_regimes] <- mean(current_reg_info[3,1:(CR_length-1)])
window_info[4,num_regimes] <- mean(current_reg_info[4,1:(CR_length-1)])
window_info[5,num_regimes] <- mean(current_reg_info[5,1:(CR_length-1)])
window_info[6,num_regimes] <- mean(current_reg_info[6,1:(CR_length-1)])
window_info[7,num_regimes] <- mean(current_reg_info[7,1:(CR_length-1)])
window_info[8,num_regimes] <- mean(current_reg_info[8,1:(CR_length-1)])
window_info[9,num_regimes] <- mean(current_reg_info[9,1:(CR_length-1)])
window_info[10,num_regimes] <- mean(current_reg_info[10,1:(CR_length-1)])
window_info[11,num_regimes] <- sd(current_reg_info[1,1:(CR_length-1)], na.rm = TRUE)
window_info[12,num_regimes] <- sd(current_reg_info[2,1:(CR_length-1)], na.rm = TRUE)
window_info[13,num_regimes] <- sd(current_reg_info[3,1:(CR_length-1)], na.rm = TRUE)
window_info[14,num_regimes] <- sd(current_reg_info[4,1:(CR_length-1)], na.rm = TRUE)
window_info[15,num_regimes] <- sd(current_reg_info[5,1:(CR_length-1)], na.rm = TRUE)
window_info[16,num_regimes] <- sd(current_reg_info[6,1:(CR_length-1)], na.rm = TRUE)
window_info[17,num_regimes] <- sd(current_reg_info[7,1:(CR_length-1)], na.rm = TRUE)
window_info[18,num_regimes] <- sd(current_reg_info[8,1:(CR_length-1)], na.rm = TRUE)
window_info[19,num_regimes] <- sd(current_reg_info[9,1:(CR_length-1)], na.rm = TRUE)
window_info[20,num_regimes] <- sd(current_reg_info[10,1:(CR_length-1)], na.rm = TRUE)
window_info[21,num_regimes] <- CR_start
window_info[22,num_regimes] <- i-1
window_info[23,num_regimes] <- cur_best_window
for (jfk in CR_start:(i-1))
{
			tw2<-xts()
			for (ii in 1:171)
			{
				tw2 <- merge(tw2, currencies[jfk:(jfk+window-1),ii])
				ii=ii+1
			}
			cor_tw2 <- cor(tw2)
			tempSVD2 <- svd(cor_tw2)
			tw2_top5 <- tempSVD2$u[,1:5]
			A <- tempSVD2$u[,1:comps_rot]
			At <- t(A)
			AtB <- At %*% W1_rot
			tempSVD3 <- svd(AtB)
			tw2_rot <- A %*% tempSVD3$u %*% t(tempSVD3$v)
			#factors_cor <= tw2_top5[,1:5] #dont want to waste time with unrotated scores
			tw2_top5 <- tw2_rot[,1:5]
			scorez <- tw2 %*% tw2_top5
			for (kk in 1:5)
			{
				factor_scores[jfk: (jfk +window-1),((jfk-1)*5+kk)] <- scorez[1:120,kk]
			}
}
# last regime has been closed out
# we have now found our regimes, and calculated the scores for the top 5 components for each window

write.csv(data_matrix, file="f2s.csv",row.names = TRUE)
Ad(get('VTI'))
