par = [0.015 0.45 0.07 -1.75 0.3 -0.06 0.45 0 0 13.63 32.16];
pt = parsecpoints(par);
evpt = evenpar(par);
plot_(pt,'r-')
plot_pt(evpt,'bo')