loadl,'fort.12',d,6,1000000
n=n_elements(d[0,*])/2
i=5
plot,d[i,0:n-1]
plot,d[i,0:n-1]-d[i,n:2*n-1]
plot,(d[i,0:n-1]-d[i,n:2*n-1])/mean(d[i,0:n-1])
l=where(d[i,0:n-1]-d[i,n:2*n-1] LT -1e-4)
print,d[0,l]
print,d[1,l]
for i=2,5 do begin & plot,d[i,0:n-1]-d[i,n:2*n-1],title=i & wait,1 & end
for i=2,5 do begin & plot,(d[i,0:n-1]-d[i,n:2*n-1])/mean(d[i,0:n-1]),title=i & wait,1 & end
