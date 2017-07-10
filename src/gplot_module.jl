module gplot_module
# This module is simply an interface with Gaston.

  using Gaston

  export gsemilogy
  export gplot
  export show_gplot
  export gfigure
  export gcloseall

  function gsemilogy(x,y,color="black")
    # set axis
    a = AxesConf(); a.axis="semilogy";   addconf(a)
    # set curve
    c = CurveConf();    c.color=color;

    # plot
    addcoords(x, y, c);
  end

  function gplot(x,y,color="black")
    # set axis
    a = AxesConf(); addconf(a)
    # set curve
    c = CurveConf();    c.color=color;

    # plot
    addcoords(x, y, c);
  end

  # call this function once you have finished plotting one figure
  show_gplot()=llplot()

  # make a new figure
  gfigure()=figure();

  # close all figures
  gcloseall()=closeall()

end
