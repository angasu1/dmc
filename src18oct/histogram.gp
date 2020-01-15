 nwidth=0.005
 bin(x,width)=width*floor(x/width)
 plot 'dynamics.dat' using (bin($2,binwidth)):(1.0) smooth freq with lines
