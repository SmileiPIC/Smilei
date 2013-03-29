set term x11 enh
i = 0

set tics out nomirr
set yr [-5:5]
set y2r [-5e-1:5e-1]
set title sprintf("%d",i*100)
bind n 'i=i+1; set title sprintf("%d",i*100) ; replot' 
bind p 'i=i-1; set title sprintf("%d",i*100) ; replot' 
 
size=0.7
width=3
set y2tics
 
plot 'fex' i i u ($0*(2.0*pi)/500):1 w l lw width tit 'E_x', \
'fey' i i u ($0*(2.0*pi)/500):1 w l lw width tit 'E_y', \
'fez' i i u ($0*(2.0*pi)/500):1 w l lw width tit 'E_z', \
'fjx' i i u ($0*(2.0*pi)/500):1 axes x1y2 w p pt 7 ps size tit 'J_x', \
'fjy' i i u ($0*(2.0*pi)/500):1 axes x1y2 w p pt 7 ps size tit 'J_y', \
'fjz' i i u ($0*(2.0*pi)/500):1 axes x1y2 w p pt 7 ps size tit 'J_z', \
'rho' i i u ($0*(2.0*pi)/500):1 axes x1y2 w p pt 7 ps size tit 'rho', \
'rho_old' i i u ($0*(2.0*pi)/500):1 axes x1y2 w p pt 7 ps size  tit 'rho_{old}'
