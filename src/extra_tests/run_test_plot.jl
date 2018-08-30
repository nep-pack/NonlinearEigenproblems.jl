cd("..")

workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory

using gplot_module


t = -2:0.01:2
st = sin(10pi*t)
ct = cos(10pi*t)
et = exp(abs(t/10))
xx = range(1, stop = 10, length = 50)
yy = 10 .^ range(-10, stop = 1, length = 50)


gcloseall()

gsemilogy(xx,yy)
show_gplot()

gfigure()
gplot(t,st,"red")
gplot(t,ct,"blue")
gplot(t,et,"green")
show_gplot()
