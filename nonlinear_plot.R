
x=seq(from=0,to=14,by=0.1)
y=4-1/(1+exp(-1.2*(x-7)))
plot(x,y,type="l")

ggplot() + geom_line(aes(x=x,y=y)) + geom_text(aes(label="M",x=-.5,y=4)) +
  geom_segment(aes(x=0,xend=0,y=3,yend=4),linetype=2) + 
  geom_segment(aes(x=7,xend=7,y=3,yend=3.5),linetype=2) + 
  geom_text(aes(label="K",x=8,y=3.5)) +
  geom_text(aes(label="L",x=.5,y=3.5)) +
  geom_text(aes(label="B",x=7,y=2.95))  +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) + theme_bw()
