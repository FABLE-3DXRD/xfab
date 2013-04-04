#!/usr/bin/env python
#Bugfix to enable ssh X11 forwarding of plots on mac
import matplotlib
matplotlib.use('GTK')
#end of bugfix
from optparse import OptionParser
import sys
import numpy as n
import pylab as pl
from xfab import symmetry,tools
#import matplotlib.axes3d as p3
import mpl_toolkits.mplot3d as p3
from ImageD11 import columnfile as ic
from matplotlib import cm,colors


if __name__=='__main__':


  def get_options():
    parser = OptionParser()
    parser.add_option("-i", "--input", action="store",
                      dest="input", type="string",
                      help="stem of input .gff")
    parser.add_option("-s", "--scale", action="store",
                      dest="scale", type="float",
                      default = None, #Do not scale grain sizes
                      help="Scale factor for grain sizes, default: 1")
    parser.add_option("-d", "--dim", action="store",
                      dest="dim", type="float",
                      default = 1, #plot -1<x,y,z<1
                      help="Plot size in mm, default: 1")
    parser.add_option("-t", "--text", action="store",
                      dest="text", type="string",
                      default = None, #do not give grain numbers on plot
                      help="Give grain numbers on plot (yes/no), default: None")
    parser.add_option("-c", "--sym", action="store",
                      dest="sym", type="string",
                      default = None, #no symmetry information
                      help="Colour according to symmetry information, options: \
                      None (default, random colours) \
                      standard (allowing all orientations in cubic space) \
                      cubic (according to cubic inverse pole figure) \
                      hexagonal (according to hexagonal inverse pole figure) \
                      s11/s22/s33/s23/s13/s12 (stress tensor components) \
                      e11/e22/e33/e23/e13/e12 (strain tensor components)")
    parser.add_option("-m", "--min", action="store",
                      dest="minimum", type="int",
                      default=None, #Use parameter-specific defaults
                      help="Give minimum for colour bar      \
                      Must be integer, use default: None to seen unit")
    parser.add_option("-M", "--max", action="store",
                      dest="maximum", type="int",
                      default=None, #Use parameter-specific defaults
                      help="Give maximum for colour bar      \
                      Must be integer, use default: None to seen unit")
    parser.add_option("-p", "--step", action="store",
                      dest="step", type="int",
                      default=None, #Use parameter-specific defaults
                      help="Give step size for colour bar      \
                      Must be integer, use default: None to seen unit")
    parser.add_option("-v", "--view", action="store",
                      dest="view", type="string",
                      default=None, help="View structure along the x/y/z axis")
    parser.add_option("-I", "--input2", action="store",
                      dest="input2", type="string",
                      default= None, help="stem of second input .gff")
    parser.add_option("-S", "--scale2", action="store",
                      dest="scale2", type="float",
                      default = None, #Do not scale grain sizes
                      help="Scale factor for grain sizes of file2, default: 1")
    parser.add_option("-T", "--text2", action="store",
                      dest="text2", type="string",
                      default = None, #do not give grain numbers on plot
                      help="Give grain numbers on plot for file2 (yes/no), default: None")
    parser.add_option("-C", "--sym2", action="store",
                      dest="sym2", type="string",
                      default = None, #no symmetry information
                      help="Colour according to symmetry information of file2         (see options for -c)")



    options , args = parser.parse_args()

    do_exit = False

    if options.input == None:
      print "\nNo stem of input .gff supplied [-i input]\n"
      do_exit = True
    if do_exit:
        parser.print_help()
        sys.exit()
    return options

    
  def read(file,sym,minimum,maximum,step):
    data = ic.columnfile('%s.gff' %file)
    #grain sizes
    if "grainno" not in data.titles:
        data.addcolumn(data.grain_id,'grainno')
    if "grainvolume" in data.titles:
        scale = max(data.grainvolume**0.3333)
        data.addcolumn(data.grainvolume**0.3333/scale*100,'size')
    else:
        data.addcolumn(100.*n.ones(data.nrows),'size')
    rodx = []
    rody = []
    rodz = []
    if "rodx" not in data.titles:
        for i in range(data.nrows):
            U = n.array([[data.U11[i],data.U12[i],data.U13[i]],
                         [data.U21[i],data.U22[i],data.U23[i]],
                         [data.U31[i],data.U32[i],data.U33[i]]])
            rod = tools.u_to_rod(U)
            rodx.append(rod[0])
            rody.append(rod[1])
            rodz.append(rod[2])
        data.addcolumn(n.array(rodx),'rodx')
        data.addcolumn(n.array(rody),'rody')
        data.addcolumn(n.array(rodz),'rodz')
    if sym == "standard":
        #grain colours, so that all orientations in Cubic space are allowed
        maxhx=62.8*3.14/180 
        minhx=-62.8*3.14/180
        maxhy=62.8*3.14/180 
        minhy=-62.8*3.14/180
        maxhz=62.8*3.14/180 
        minhz=-62.8*3.14/180; 
        rr = data.rodx**2+data.rody**2+data.rodz**2
        rr = n.sqrt(rr)
        theta=2*n.arctan(rr)
        r1=data.rodx/rr 
        r2=data.rody/rr 
        r3=data.rodz/rr
        # normalise colours
        red =   (r1*theta-minhx)/(maxhx-minhx)
        green = (r2*theta-minhy)/(maxhy-minhy)
        blue =  (r3*theta-minhz)/(maxhz-minhz)
    elif sym == "cubic":
        # Fill cubic stereographic triangle
        red = []
        green = []
        blue = []
        fig = pl.figure(10,frameon=False)
        ax = pl.Axes(fig,[.2,.2,.7,.7])
        ax.set_axis_off()
        fig.add_axes(ax)
        #plot triangle    
        xa = n.zeros((21))
        ya = n.zeros((21))
        for i in range(21):
            ua = n.array([i/20., 1., 1.])
            UA = n.linalg.norm(ua)
            za = ua[2]+UA
            xa[i] = ua[1]/za
            ya[i] = ua[0]/za
        pl.plot(xa,ya,'black') # Curved edge
        pl.plot([0,xa[0]],[0,0.0001],'black',linewidth=2) #lower line 
        pl.plot([xa[20],0],[ya[20],0],'black') # upper line
        pl.text(-0.01,-0.02,'[100]')
        pl.text(xa[0]-0.01,-0.02,'[110]')
        pl.text(xa[-1]-0.01,ya[-1]+0.005,'[111]')
        # Grains
        for i in range(data.nrows):
            U = tools.rod_to_u([data.rodx[i],data.rody[i],data.rodz[i]])
            axis = abs(U[2,:])
            colour = n.zeros((3))
            for j in range(3):
                for k in range(j+1,3):
                    if (axis[j]>axis[k]):
                        colour[0]=axis[j]
                        axis[j]=axis[k]
                        axis[k]=colour[0]                 
            rr=n.sqrt(axis[0]*axis[0]/((axis[2]+1))/((axis[2]+1))+(axis[1]/(axis[2]+1)+1)*(axis[1]/(axis[2]+1)+1))
            if axis[1]==0: 
                beta=0
            else:
                beta=n.arctan(axis[0]/axis[1])
            colour[0]=((n.sqrt(2.0)-rr)/(n.sqrt(2.0)-1))**.5;
            colour[1]=((1-4*beta/n.pi)*((rr-1)/(n.sqrt(2.0)-1)))**.5;
            colour[2]=(4*beta/n.pi*((rr-1)/(n.sqrt(2.0)-1)))**.5;
            red.append(colour[0])
            green.append(colour[1])
            blue.append(colour[2])
            X = axis[1]/(1+axis[2])
            Y = axis[0]/(1+axis[2])            
            pl.plot(X,Y,'o',color=colour)
        pl.xlim()
    elif sym == "hexagonal":
        # Fill hexagonal stereographic triangle
        A = n.array([[0,1,0],[0,-n.sqrt(3),1],[1,0,0]])
        a0 = 1./n.sqrt(3.)
        a1 = 1.
        a2 = 1.
        red = []
        green = []
        blue = []
        fig = pl.figure(10,frameon=False)
        ax = pl.Axes(fig,[0.2,.2,0.6,0.6])
        ax.set_axis_off()
        fig.add_axes(ax)
        #plot triangle    
        ya = n.array(range(51))/100.
        xa = n.sqrt(1-ya**2)
        pl.plot(xa,ya,'black') # Curved edge
        pl.plot([0,xa[0]],[0,0.0001],'black',linewidth=2) #lower line 
        pl.plot([xa[-1],0],[ya[-1],0],'black') # upper line
        pl.text(-0.01,-0.02,'[001]')
        pl.text(xa[0]-0.03,-0.02,'[010]')
        pl.text(xa[-1]-0.03,ya[-1]+0.005,'[120]')
        # Grains
        r = symmetry.rotations(6)
        for i in range(data.nrows):
            U = tools.rod_to_u([data.rodx[i],data.rody[i],data.rodz[i]])
            square = 1
            angle = 0
            frac = 1./n.sqrt(3.)
            for k in range(len(r)):
                g = n.dot(U,r[k])
                a = n.arccos((n.trace(g)-1)/2)                
                if g[2,2] > 0:
                    uvw = g[2,:]
                else:
                    uvw = -g[2,:]
                x = uvw[0]/(1+uvw[2])
                y = uvw[1]/(1+uvw[2])
                f = y/x
                s = x*x+y*y
                if f<=frac and s<=square and x>=0 and y>=0:
                    angle = a
                    frac = f
                    square = s
                    UVW = uvw
                    X = x
                    Y = y
            color = n.dot(n.transpose(A),n.transpose(UVW))**0.7
#            print color
            color[0] = color[0]/a0
            color[1] = color[1]/a1
            color[2] = color[2]/a2
            mx = max(color)
            color = color/mx
            red.append(color[0])
            green.append(color[1])
            blue.append(color[2])            
            pl.plot(X,Y,'o',color=color)
        pl.xlim()
    elif sym == "e11":
        try:
            colourbar(minimum,maximum,step,'e-3')
            norm = colors.normalize(minimum*1e-3,maximum*1e-3)
        except:
            colourbar(-1,1,1,'e-3')
            norm = colors.normalize(-1e-3,1e-3)
        color = cm.jet(norm(data.eps11_s))
        red = color[:,0]
        green = color[:,1]
        blue = color[:,2]
        print "transverse strain:", n.sum(data.grainvolume*data.eps11_s)/n.sum(data.grainvolume)
    elif sym == "e22":
        try:
            colourbar(minimum,maximum,step,'e-3')
            norm = colors.normalize(minimum*1e-3,maximum*1e-3)
        except:
            colourbar(-1,1,1,'e-3')
            norm = colors.normalize(-1e-3,1e-3)
        color = cm.jet(norm(data.eps22_s))
        red = color[:,0]
        green = color[:,1]
        blue = color[:,2]
        print "transverse strain:", n.sum(data.grainvolume*data.eps22_s)/n.sum(data.grainvolume)
    elif sym == "e33":
        try:
            colourbar(minimum,maximum,step,'e-3')
            norm = colors.normalize(minimum*1e-3,maximum*1e-3)
        except:
            colourbar(-1,1,1,'e-3')
            norm = colors.normalize(-1e-3,1e-3)
        color = cm.jet(norm(data.eps33_s))
        red = color[:,0]
        green = color[:,1]
        blue = color[:,2]
        print "axial strain:", n.sum(data.grainvolume*data.eps33_s)/n.sum(data.grainvolume)
    elif sym == "e12":
        try:
            colourbar(minimum,maximum,step,'e-3')
            norm = colors.normalize(minimum*1e-3,maximum*1e-3)
        except:
            colourbar(-1,1,1,'e-3')
            norm = colors.normalize(-1e-3,1e-3)
        color = cm.jet(norm(data.eps12_s))
        red = color[:,0]
        green = color[:,1]
        blue = color[:,2]
        print "shear strain:", n.sum(data.grainvolume*data.eps12_s)/n.sum(data.grainvolume)
    elif sym == "e13":
        try:
            colourbar(minimum,maximum,step,'e-3')
            norm = colors.normalize(minimum*1e-3,maximum*1e-3)
        except:
            colourbar(-1,1,1,'e-3')
            norm = colors.normalize(-1e-3,1e-3)
        color = cm.jet(norm(data.eps13_s))
        red = color[:,0]
        green = color[:,1]
        blue = color[:,2]
        print "shear strain:", n.sum(data.grainvolume*data.eps13_s)/n.sum(data.grainvolume)
    elif sym == "e23":
        try:
            colourbar(minimum,maximum,step,'e-3')
            norm = colors.normalize(minimum*1e-3,maximum*1e-3)
        except:
            colourbar(-1,1,1,'e-3')
            norm = colors.normalize(-1e-3,1e-3)
        color = cm.jet(norm(data.eps23_s))
        red = color[:,0]
        green = color[:,1]
        blue = color[:,2]
        print "shear strain:", n.sum(data.grainvolume*data.eps23_s)/n.sum(data.grainvolume)
    elif sym == "s33":
        try:
            colourbar(minimum,maximum,step,'MPa')
            norm = colors.normalize(minimum,maximum)
        except:
            colourbar(-50,150,50,'MPa')
            norm = colors.normalize(-50,150)
        color = cm.jet(norm(data.sig33_s))
        red = color[:,0]
        green = color[:,1]
        blue = color[:,2]
        print "axial stress:", n.sum(data.grainvolume*data.sig33_s)/n.sum(data.grainvolume), "MPa" 
    elif sym == "s11":
        try:
            colourbar(minimum,maximum,step,'MPa')
            norm = colors.normalize(minimum,maximum)
        except:
            colourbar(-50,150,50,'MPa')
            norm = colors.normalize(-50,150)
        color = cm.jet(norm(data.sig11_s))
        red = color[:,0]
        green = color[:,1]
        blue = color[:,2]
        print "transverse s11 stress:", n.sum(data.grainvolume*data.sig11_s)/n.sum(data.grainvolume), "MPa" 
    elif sym == "s22":
        try:
            colourbar(minimum,maximum,step,'MPa')
            norm = colors.normalize(minimum,maximum)
        except:
            colourbar(-50,150,50,'MPa')
            norm = colors.normalize(-50,150)
        color = cm.jet(norm(data.sig22_s))
        red = color[:,0]
        green = color[:,1]
        blue = color[:,2]
        print "transverse s22 stress:", n.sum(data.grainvolume*data.sig22_s)/n.sum(data.grainvolume), "MPa" 
    elif sym == "s12":
        try:
            colourbar(minimum,maximum,step,'MPa')
            norm = colors.normalize(minimum,maximum)
        except:
            colourbar(-50,50,50,'MPa')
            norm = colors.normalize(-50,50)
        color = cm.jet(norm(data.sig12_s))
        red = color[:,0]
        green = color[:,1]
        blue = color[:,2]
        print "shear s12 stress:", n.sum(data.grainvolume*data.sig12_s)/n.sum(data.grainvolume), "MPa" 
    elif sym == "s13":
        try:
            colourbar(minimum,maximum,step,'MPa')
            norm = colors.normalize(minimum,maximum)
        except:
            colourbar(-50,50,50,'MPa')
            norm = colors.normalize(-50,50)
        color = cm.jet(norm(data.sig13_s))
        red = color[:,0]
        green = color[:,1]
        blue = color[:,2]
        print "shear s13 stress:", n.sum(data.grainvolume*data.sig13_s)/n.sum(data.grainvolume), "MPa" 
    elif sym == "s23":
        try:
            colourbar(minimum,maximum,step,'MPa')
            norm = colors.normalize(minimum,maximum)
        except:
            colourbar(-50,50,50,'MPa')
            norm = colors.normalize(-50,50)
        color = cm.jet(norm(data.sig23_s))
        red = color[:,0]
        green = color[:,1]
        blue = color[:,2]
        print "shear s23 stress:", n.sum(data.grainvolume*data.sig23_s)/n.sum(data.grainvolume), "MPa" 
    elif sym == "latt_rot":
        norm = colors.normalize(0,0.5)
        color = cm.jet(norm(data.latt_rot))
        red = color[:,0]
        green = color[:,1]
        blue = color[:,2]
        colourbar(0.0,0.5,0.1,'deg')
    elif sym == "tz":
        norm = colors.normalize(-0.1,0.1)
        color = cm.jet(norm(data.tz))
        red = color[:,0]
        green = color[:,1]
        blue = color[:,2]
    elif sym == "vol":
        norm = colors.normalize(0,10)
        color = cm.jet(norm(data.grainvolume/data.d_grainvolume))
        red = color[:,0]
        green = color[:,1]
        blue = color[:,2]
    elif sym == "tth":
        norm = colors.normalize(0.007,0.009)
        color = cm.jet(norm(data.sig_tth/data.grainvolume**.2))
        red = color[:,0]
        green = color[:,1]
        blue = color[:,2]
        print min(data.sig_tth/data.grainvolume**.2),max(data.sig_tth/data.grainvolume**.2)
#        pl.figure(8)
#        pl.plot(n.log(data.grainvolume),n.log(data.sig_tth),'.')
    elif sym == "eta":
        norm = colors.normalize(0.08,0.15)
        color = cm.jet(norm(data.sig_eta/data.grainvolume**.2))
        red = color[:,0]
        green = color[:,1]
        blue = color[:,2]
        print min(data.sig_eta/data.grainvolume**.2),max(data.sig_eta/data.grainvolume**.2)
#        pl.figure(8)
#        pl.plot(n.log(data.grainvolume),n.log(data.sig_eta),'.')
    else:
        n.random.seed(0)
        red = n.random.rand(data.nrows)
        n.random.seed(1)
        green = n.random.rand(data.nrows)
        n.random.seed(2)
        blue = n.random.rand(data.nrows)
    data.addcolumn(red,'red')
    data.addcolumn(green,'green')
    data.addcolumn(blue,'blue')
    return (data)
 
 
  def scale(data,scale):
    """
    scale grain sizes
    """
    if scale != None:
        data.addcolumn(data.size*scale,'size')
            
  
  def colourbar(min,max,step,legend):
    """
    add colourbar to plot with min, max and legend
    """
    fig2 = pl.figure(0)
    ax = fig2.add_subplot(111,visible=False)
    data = n.clip(n.random.randn(250,250),-1,1)
    cax = ax.imshow(data,interpolation='nearest',cmap=cm.jet)
    numticks = range(min,max+step,step)
    textticks = []
    num2ticks = []
    for i in numticks:
        textticks.append('%s' %i)
        num2ticks.append(2.*(i-(min+max)/2.)/(max-min))
#    print num2ticks
#    print textticks
    cbar = fig2.colorbar(cax,ticks=num2ticks,orientation='horizontal')
    cbar.ax.set_xticklabels(textticks)
    cbar.ax.set_title('[%s]' %legend,va='baseline')
  
  
  def plot3Dgrains(data,dim,view=None,label=None,data2=None,label2=None):
    """
    Plots xyz of grains in 3D
    Jette Oddershede March 4th 2009
    """
    if view == "z":
        elev = 90
        azim = 270
    elif view == "y":
        elev = 0
        azim = 270
    else:    
        elev = 0
        azim = 0
    fig = pl.figure()
    ax = p3.Axes3D(fig)
    #plot data
    color = []
    for i in range(data.nrows):
        color.append([data.red[i],data.green[i],data.blue[i]])
    cax=ax.scatter3D(data.x,data.y,data.z,s=data.size,c=color)
    if data2 != None:
        color2 = []
        for i in range(data2.nrows):
            color2.append([data2.red[i],data2.green[i],data2.blue[i]])
        ax.scatter3D(data2.x,data2.y,data2.z,s=data2.size,c=color)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_xlim3d(-dim,dim)
    ax.set_ylim3d(-dim,dim)
    ax.set_zlim3d(-dim,dim)
    ax.view_init(elev, azim) 
    if label != None:
        for i in range(data.nrows):
            ax.text(data.x[i],data.y[i],data.z[i],'%i' %(data.grainno[i]))
    if data2 != None:
        if label2 != None:
           for i in range(data2.nrows):
                ax.text(data2.x[i],data2.y[i],data2.z[i],'%i' %(data2.grainno[i]))
    pl.title('%s    %s'%(options.input,options.sym))
    pl.show()
    
  options = get_options()
  (data) = read(options.input,options.sym,options.minimum,options.maximum,options.step)
  if options.input2 != None:
    (data2) = read(options.input2,options.sym2)  
    scale(data2,options.scale2)
  scale(data,options.scale)
  if options.text == "no" or options.text == "NO" or options.text == "No" or options.text == "n" or options.text == "N":
    options.text = None
  if options.input2 == None:
    plot3Dgrains(data,options.dim,options.view,label=options.text)
  else:
    plot3Dgrains(data,options.dim,options.view,label=options.text,data2=data2,label2=options.text2)

    

    


    

    
