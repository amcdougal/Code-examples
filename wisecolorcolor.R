#This program plots the W1-W2 vs W2-W3 colors for the allsky WISE-FIRST
# catalog.  Plots them as contours.  Only use sources with SNR>5 in the
# relevant band. 
# Creates Figures 7 & 8 in Truebenbach & Darling, 2017, MNRAS, 468, 196 
#By: Alex Truebenbach
#Last Editted: Jan 2016

#Needed for kde2d
library(MASS)


#Read in file
allsky=read.csv("../final_match_allskyV.no2mass_FINAL.short.tbl")

#Select those with SNR>5 in the first 3 bands and exclude extended sources
highw1w2w3=allsky$snr1>=5 & allsky$snr2>=5 & allsky$snr3>=5 & allsky$ex<=0

#Remove those with NA snr values
newsnr1=allsky$snr1[highw1w2w3]
newsnr2=allsky$snr2[highw1w2w3]
newsnr3=allsky$snr3[highw1w2w3]
no.na=c(which(complete.cases(newsnr1)),which(complete.cases(newsnr2)),which(complete.cases(newsnr3)))

no.na=unique(no.na)

w1=allsky$W1mag[highw1w2w3]
w2=allsky$W2mag[highw1w2w3]
w3=allsky$W3mag[highw1w2w3]
ra=allsky$search_ra[highw1w2w3]
dec=allsky$search_dec[highw1w2w3]
ra=ra[no.na]
dec=dec[no.na]

#Find w1-w2 and w2-w3 colors for objects detected in the first 3 bands
w1w2=w1[no.na]-w2[no.na]
w2w3=w2[no.na]-w3[no.na]

print("Eligible for Mateos wedge:")
print(length(w1w2))
print("")

#Calculate the 2d kernel density
z=kde2d(w2w3,w1w2,n=100)

#Plotting
pdf(file="../../paper_graphs/colorcolor.pdf",height=10,width=10)
#Plot the outlying points
smoothScatter(w2w3,w1w2,nrpoints=320,nbin=200,cex=.8,pch=19,colramp=colorRampPalette(c("white","white")), xlim=c(1.5,6.5),ylim=c(0,3),xlab="[4.6]-[12]",ylab="[3.4]-[4.6]",cex.lab=1.5,cex.axis=1.5)
#Plot the contour lines
contour(z,drawlabels=F,nlevels=5,add=T,lwd=2)

#Add Stern+2012 w1-w2>0.8 line
abline(h=0.8,lty="dashed",col="red",lwd=2)
text(6,.9,"Stern+12",col="red",cex=1.5)

#Add Mateos+ 2012 Wedge
#Top line:
segments(1.958,1.412,7,3.001,lty="dotted",col="blue",lwd=2.5)
#Bottom line:
segments(2.25,0.487,7,1.983,lty="dotted",col="blue",lwd=2.5)
#Left line:
segments(1.958,1.412,2.25,0.487,lty="dotted",col="blue",lwd=2.5)
text(6,2.4,"Mateos+12",col="blue",cex=1.5)

#Add in z evolution of Arp220 and Pope SMG template
smg=read.csv("/duanestorage/home/student/altr4657/Documents/research/WISEcolors/WISEcolors_smg_pope08.tbl")
arp=read.csv("/duanestorage/home/student/altr4657/Documents/research/WISEcolors/WISEcolors_Arp220.tbl")
#convert to Vega mags
smgW1=(log10(smg$W1[1:51]/3631)*(-5/2))-2.7
smgW2=log10(smg$W2[1:51]/3631)*(-5/2)-3.34
smgW3=log10(smg$W3[1:51]/3631)*(-5/2)-5.17
arpW1=(log10(arp$W1[1:51]/3631)*-2.5)-2.7
arpW2=log10(arp$W2[1:51]/3631)*(-5/2)-3.34
arpW3=log10(arp$W3[1:51]/3631)*(-5/2)-5.17

#Plot
x=smgW2-smgW3
y=smgW1-smgW2

#Create color gradient
colors=colorRampPalette(c("palegreen","darkgreen"))
arrows(x[1:(length(x)-1)],y[1:(length(x)-1)],x[2:length(x)],y[2:length(x)],length=.0001,col=colors(50),lwd=3)

#Add in markers for integer redshifts
points(c(x[11],x[21],x[31],x[41]),c(y[11],y[21],y[31],y[41]),col="black", bg=c("palegreen2","springgreen","seagreen3","springgreen4"),cex=3.5,pch=23)
 points(c(x[11],x[21],x[31],x[41]),c(y[11],y[21],y[31],y[41]),col="black",pch=c("1","2","3","4"),cex=2)

#Find redshifts where smg is selected by wedge/stern cut
print("Selected Redshifts for SMG:")
print("Mateos Wedge:")
d <- c()
for (index in 1:length(x)) {
    a=0.315*x[index]-.222
    b=(y[index]-7.624)/-3.172
    if ((y[index]>a) * (x[index]>b)) {
       d<-c(d,smg$redshift[index])}}
print(d)

print("Stern Line:")
print(smg$redshift[y>0.8])

x=arpW2-arpW3
y=arpW1-arpW2
colors=colorRampPalette(c("orchid","darkorchid4"))
arrows(x[1:(length(x)-1)],y[1:(length(x)-1)],x[2:length(x)],y[2:length(x)],length=.0001,col=colors(50),lwd=3)


#Add in markers for integer redshifts
points(c(x[11],x[21],x[31],x[41]),c(y[11],y[21],y[31],y[41]),col="black", bg=c("orchid","mediumorchid","darkviolet","darkorchid4"),cex=3.5,pch=23)
 points(c(x[11],x[21],x[31],x[41]),c(y[11],y[21],y[31],y[41]),col="black",cex=2,pch=c("1","2","3","4"))

#Find redshifts where smg is selected by wedge/stern cut
print("Selected Redshifts for Arp220:")
print("Mateos Wedge:")
d <- c()
for (index in 1:length(x)) {
    a=0.315*x[index]-.222
    b=(y[index]-7.624)/-3.172
    if ((y[index]>a) * (x[index]>b)) {
       d<-c(d,smg$redshift[index])}}
print(d)

print("Stern Line:")
print(smg$redshift[y>0.8])

#Legend
legend(1.5,3,c("SMG","Arp 220"),col=c("seagreen3","mediumorchid"),lwd=2,cex=1.5)

dev.off()



##Calc number in mateos wedge
#1 is to left of wedge.  Besides that, calc number below top line - number
#below bottom line.  -1 for 1 to the left of the wedge.
print(" ")
print("Numbers for calculating total selected:")
print("Mateos wedge:")
n=0
for (index in 1:length(w1w2)) {
    if ((w1w2[index]-1.412)<(0.315*(w2w3[index]-1.958))) {
       n=n+1}
       }

print("below top line")
print(n)

n=0
for (index in 1:length(w1w2)) {
    if ((w1w2[index]-0.487)<(0.315*(w2w3[index]-2.25))) {
       n=n+1}
       }

print("below bottom line")
print(n)


###Also plot W1-W2 by itself.  Similar to Stern+2012
highw1w2=allsky$snr1>5 & allsky$snr2>5 & allsky$ex<=0
#Remove those with NA snr values
newsnr1=allsky$snr1[highw1w2]
newsnr2=allsky$snr2[highw1w2]
no.na=c(which(complete.cases(newsnr1)),which(complete.cases(newsnr2)))

no.na=unique(no.na)

w1=allsky$W1mag[highw1w2]
w2=allsky$W2mag[highw1w2]

w1w2=w1[no.na]-w2[no.na]
w1=w1[no.na]

print("Eligible for stern:")
print(length(w1w2))

#Plot the outlying points
pdf(file="../../paper_graphs/w1w2color.pdf",height=10,width=10)
smoothScatter(w1,w1w2,nrpoints=4500,nbin=200,cex=.5,pch=19,colramp=colorRampPalette(c("white","white")),xlab="[3.4]",ylab="[3.4]-[4.6]",ylim=c(-.5,3),xlim=c(13,19),cex.lab=1.5,cex.axis=1.5)
#Calculate the 2d kernel density
z=kde2d(w1,w1w2,n=50)
#Plot the contour lines
contour(z,drawlabels=T,add=T,lwd=2)

#Add Stern+2012 w1-w2>0.8 line
abline(h=0.8,lty="dashed",col="red",lwd=2)
#text(0,1,"Stern+12",col="red",cex=1.5)


#Add Arp220 and SMG lines
x=smgW1
y=smgW1-smgW2
#Create color gradient
colors=colorRampPalette(c("springgreen","darkgreen"))
arrows(x[4:(length(x)-11)],y[4:(length(x)-11)],x[5:(length(x)-10)],y[5:(length(x)-10)],length=.0001,col=colors(49),lwd=3)
colors=colorRampPalette(c("palegreen","darkgreen"))
arrows(x[1:(length(x)-1)],y[1:(length(x)-1)],x[2:length(x)],y[2:length(x)],length=.0001,col=colors(50),lwd=3)

#Add in markers for integer redshifts
points(c(x[11],x[21],x[31],x[41]),c(y[11],y[21],y[31],y[41]),col="black", bg=c("palegreen2","springgreen","seagreen3","springgreen4"),cex=3.5,pch=23)
 points(c(x[11],x[21],x[31],x[41]),c(y[11],y[21],y[31],y[41]),col="black",pch=c("1","2","3","4"),cex=2)

#smg from z=.4-4

x=arpW1
y=arpW1-arpW2
print(x)
colors=colorRampPalette(c("orchid","darkviolet"))
arrows(x[2:(length(x)-35)],y[2:(length(x)-35)],x[3:(length(x)-34)],y[3:(length(x)-34)],length=.0001,col=colors(24),lwd=3)
#Add in markers for integer redshifts
points(c(x[11]),c(y[11]),col="black", bg=c("orchid"),cex=3.5,pch=23)
 points(c(x[11]),c(y[11]),col="black",pch=c("1"),cex=2)

#Arp220 z from 0.1-1.6


#Legend
legend(13,3,c("SMG","Arp 220"),col=c("seagreen3","mediumorchid"),lwd=2,cex=1.5)
ra=allsky$search_ra[highw1w2]
dec=allsky$search_dec[highw1w2]
ra=ra[no.na]
dec=dec[no.na]

#WISE 95% compelete for mags brighter than 17.1 in W1. Varies across sky
# depending on survey coverage, backround emission, and source density.
#Add this line
abline(v=17.1,col="blue",lwd=2.5,lty=3)

#Calc number w/ W1 and W2 detections:
print("W1-W2 detections:")
print(length(w1w2))
print(length(w1[w1w2>0.8]))

dev.off()
