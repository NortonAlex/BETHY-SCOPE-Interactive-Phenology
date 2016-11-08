PRO line, dline=dline

infile='line.set'
IF KEYWORD_SET(dline) THEN infile='dline.set'
loadl,infile,d,2,100000
set_plot,'ps'
psfile='line.ps'
IF KEYWORD_SET(dline) THEN psfile='dline.ps'
device,file=psfile
ytitle='Delta CF'
IF KEYWORD_SET(dline) THEN ytitle='Delta gradient CF'
plot,d[0,*],d[1,*],xtitle='norm of difference of control vector to base point',ytitle=ytitle
device,/close

END

