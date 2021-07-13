install.packages("hexSticker")
library(hexSticker)

png("inst/figures/zalpha_image.png",res=100,bg="transparent",width=15,height=10,units = "cm")
par(mar=c(1,1,1,1))
plot_x<-seq(0,15,0.01)
plot_norm<-seq(-4,4,length.out=length(plot_x))
plot_line_1<-sin(plot_x*pi)+10*dnorm(plot_norm)
plot_line_2<--sin(plot_x*pi)+10*dnorm(plot_norm)

plot(plot_x,plot_line_1,type="n",ylim=c(-1,6),lwd=6,xaxt="n",yaxt="n",ylab="",xlab="",frame.plot = FALSE)
for(i in 0:14){
  segments(plot_x[34+i*100], plot_line_1[34+i*100], y1 = plot_line_2[34+i*100], lwd=3)
  segments(plot_x[66+i*100], plot_line_1[66+i*100], y1 = plot_line_2[66+i*100], lwd=3)
}
lines(plot_x[1:51],plot_line_1[1:51],col="black",lwd=6)
lines(plot_x[1:51],plot_line_2[1:51],col="blue",lwd=6)
for(i in 0:6){
  lines(plot_x[c(52:151)+i*200],plot_line_1[c(52:151)+i*200],col="black",lwd=6)
  lines(plot_x[c(52:151)+i*200],plot_line_2[c(52:151)+i*200],col="blue",lwd=6)
}
for(i in 0:6){
  lines(plot_x[c(152:251)+i*200],plot_line_2[c(152:251)+i*200],col="blue",lwd=6)
  lines(plot_x[c(152:251)+i*200],plot_line_1[c(152:251)+i*200],col="black",lwd=6)
}
lines(plot_x[1451:1501],plot_line_1[1451:1501],col="black",lwd=6)
lines(plot_x[1451:1501],plot_line_2[1451:1501],col="blue",lwd=6)

axis(side=1,lwd=2)
axis(side=2,at=seq(-1,3,1),lwd=2)
dev.off()

imgurl <- "inst/figures/zalpha_image.png"
sticker(imgurl, package="zalpha", p_size=30, s_x=1, s_y=1.4, s_width=.9, p_y=0.6,
        h_fill="red",h_color = "black",spotlight=TRUE,l_x=1,l_y=1,#l_w=4,l_h=4,l_alpha=0.5,
        filename="inst/figures/sticker.png")
