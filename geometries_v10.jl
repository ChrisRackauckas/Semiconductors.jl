# Semiconductor geometry collection
# Sam Chinnery
# 2022-03-09

################################################################################
# PRIMITIVES
################################################################################

# Create an arc centered at c, with radius r, from angle t1 to t2 (in degrees)
# This assumes that the appropriate facetregion is already set
# It also assumes that the endpoints are already in the grid as p1 and p2
function arc!(b,r,c,t1,t2,p1,p2;N=10,t_segment=nothing)

	# Convert angles to radians
	@assert t1<t2
	deg_to_rad = pi/180
	t1_rad = t1*deg_to_rad
	t2_rad = t2*deg_to_rad

	# Create points
	if t_segment==nothing
		t_range = range(t1_rad,t2_rad,length=N)
		points = [point!(b,c[1]+r*cos(t),c[2]+r*sin(t)) 
			for t in t_range[2: end-1]]
	else
		t_segment_rad = t_segment*pi/180
		t_range = t1_rad: t_segment_rad: t2_rad
		points = [point!(b,c[1]+r*cos(t),c[2]+r*sin(t))
			for t in t_range[2: end-1]]
		N = length(points)
	end

	# Create facets
	for i in 1: N-3
		facet!(b,points[i],points[i+1])
	end
	facet!(b,p1,points[1])
	facet!(b,points[end],p2)
	
end

# Create a rectangular background mesh using the specified grid builder
# The background mesh should not share a boundary with any existing region
function bgmesh!(b,b_types,r_types,b_num,r_num,x1,x2,y1,y2,max_vol;
	builder=false)

	# Background mesh
	@assert x1<x2
	@assert y1<y2
	p1 = point!(b,x1,y1)
	p2 = point!(b,x1,y2)
	p3 = point!(b,x2,y2)
	p4 = point!(b,x2,y1)

	# Exterior facets
	facetregion!(b,b_num)
	facet!(b,p1,p2)
	facet!(b,p2,p3)
	facet!(b,p3,p4)
	facet!(b,p4,p1)

	# Background mesh region
	# TODO: Find a better way to pick the regionpoint
	cellregion!(b,r_num)
	regionpoint!(b,x1*0.999,y1*0.999)
	options!(b,maxvolume=max_vol)

	# Add new boundary and region to type dictionaries
	merge!(b_types,Dict(b_num=>"background"))
	merge!(r_types,Dict(r_num=>"insulator"))

	# Return mesh
	if !builder
		return (simplexgrid(b),b_types,r_types)
	else
		return (b,b_types,r_types)
	end

end

################################################################################
# DIODES
################################################################################

# 1-D planar diode with contacts at either end
function diode1d_grid(wp,wn,ha,hb,hc)

	# Create grid refined about junction center
	grid_p = geomspace(-wp,0,ha,hb)
	grid_n = geomspace(0,wn,hb,hc)[2:end]
	grid_list = vcat(grid_p,grid_n)
	grid = simplexgrid(grid_list)

	# Create P and N-type regions
	cellmask!(grid,[-wp],[0],1)
	cellmask!(grid,[0],[wn],2)
	bfacemask!(grid,[0],[0],3)

	# Boundary and region types
	b_types = Dict(1=>"contact",2=>"contact",3=>"interior")
	r_types = Dict(1=>"semiconductor",2=>"semiconductor")

	return (grid,b_types,r_types)

end

# 1-D planar diode with contacts at either end
function diode1d(wp,wn;ha=0.05,hb=0.01,hc=0.05,na=1e3,nd=1e3,doping="abrupt")

	grid,b_types,r_types = diode1d_grid(wp,wn,ha,hb,hc)
	model = Model()

	# Material properties
	model.er = Dict(1=>11.7,2=>11.7)
	model.ni = Dict(1=>1.1e-2,2=>1.1e-2)
	model.ni_boundary = Dict(1=>1.1e-2,2=>1.1e-2,3=>1.1e-2)
	model.tn = Dict(1=>1e-10,2=>1e-10)
	model.tp = Dict(1=>1e-10,2=>1e-10)

	# User-selectable doping profile
	function doping_abrupt(x,y,reg;boundary=false)
		if !boundary
			if reg==1
				return (-na)
			else
				return (nd)
			end
		else
			if reg==1
				return (-na)
			elseif reg==2
				return (nd)
			else
				return (0)
			end
		end
	end
	function doping_lin(x,y,reg;boundary=false)
		return ( (wp*nd-wn*na+(nd+na)*x) / (wp+wn) )
	end
	function doping_smooth(x,y,reg;boundary=false)
		a = 1.0e1
		b = (nd-na)/(nd+na)
		return ( 1/2 * (nd+na) * (tanh(a*x-atanh(b))+b) )
	end
	if doping=="abrupt"
		model.doping = doping_abrupt
	elseif doping=="lin"
		model.doping = doping_lin
	elseif doping=="smooth"
		model.doping = doping_smooth
	end

	return (Device(model,grid,b_types,r_types))

end

# 2-D planar diode with contacts at either end
function diode2d_grid(wn,wp,h,hc,max_vol,refine_param,yplane;builder=false)

	b = SimplexGridBuilder(Generator=Triangulate)

	# Device dimensions
	yc = h/2
	y1 = yc-hc/2
	y2 = yc+hc/2

	# Grid points
	p1 = point!(b,-wp,0)
	p2 = point!(b,-wp,y1)
	p3 = point!(b,-wp,y2)
	p4 = point!(b,-wp,h)
	p5 = point!(b,0,h)
	p6 = point!(b,wn,h)
	p7 = point!(b,wn,y2)
	p8 = point!(b,wn,y1)
	p9 = point!(b,wn,0)
	p10 = point!(b,0,0)

	# p-side contact facet
	facetregion!(b,1)
	facet!(b,p2,p3)

	# n-side contact facet
	facetregion!(b,2)
	facet!(b,p7,p8)

	# Non-contact facets
	facetregion!(b,3)
	facet!(b,p1,p2)
	facet!(b,p3,p4)
	facet!(b,p4,p6)
	facet!(b,p6,p7)
	facet!(b,p8,p9)
	facet!(b,p9,p1)

	# Interior facets
	facetregion!(b,4)
	facet!(b,p10,p5)

	# Boundary and region types
	b_types = Dict(1=>"contact",2=>"contact",3=>"exterior",4=>"interior")
	r_types = Dict(1=>"semiconductor",2=>"semiconductor")

	# Optional cut plane along y axis
	if yplane!=nothing

		# Create facet for the cut plane
		p11 = point!(b,-wp,yplane)
		p12 = point!(b,wn,yplane)
		facetregion!(b,5)
		facet!(b,p11,p12)
		merge!(b_types,Dict(5=>"interior"))

		# Divide into P and N regions, including cut plane
		cellregion!(b,1)
		regionpoint!(b,-wp/2,yplane/2)
		regionpoint!(b,-wp/2,(h+yplane)/2)
		cellregion!(b,2)
		regionpoint!(b,wn/2,yplane/2)
		regionpoint!(b,wn/2,(h+yplane)/2)

	else

		# Divide into P and N regions without cut plane
		cellregion!(b,1)
		regionpoint!(b,-wp/2,h/2)
		cellregion!(b,2)
		regionpoint!(b,wn/2,h/2)

	end

	# Refine about middle of junction
	function unsuitable(x1,y1,x2,y2,x3,y3,area)
		xcenter = (x1+x2+x3)/3
		ycenter = (y1+y2+y3)/3
		dist = abs(xcenter)
		if (area^2 > refine_param*dist) && (-wp<=xcenter<=wn) &&
			(0<=ycenter<=h)
			return 1
		else
			return 0
		end
	end

	# Create grid
	options!(b,unsuitable=unsuitable)
	options!(b,maxvolume=max_vol)
	if !builder
		return (simplexgrid(b),b_types,r_types)
	else
		return (b,b_types,r_types)
	end

end

# 2-D planar diode with contacts at either end, rectangular grid
function diode2d_grid_rect(wn,wp,h,hc,dx1,dx2,dx3,ny)

	# Create grid arrays
	x1 = geomspace(-wp,0.0,dx1,dx2)
	x2 = geomspace(0.0,wn,dx2,dx3)
	x = vcat(x1,x2[2: end])
	y = range(0.0,h,length=ny)
	grid = simplexgrid(x,y)

	# Contact regions
	# NOTE: Facets have to be drawn in increasing x and y order to work
	bfacemask!(grid,[-wp,0.5*(h-hc)],[-wp,0.5*(h+hc)],1)
	bfacemask!(grid,[wn,0.5*(h-hc)],[wn,0.5*(h+hc)],2)

	# Exterior non-contact facets
	bfacemask!(grid,[-wp,0.0],[-wp,0.5*(h-hc)],3)
	bfacemask!(grid,[-wp,0.5*(h+hc)],[-wp,h],3)
	bfacemask!(grid,[-wp,h],[0.0,h],3)
	bfacemask!(grid,[0.0,h],[wn,h],3)
	bfacemask!(grid,[wn,0.5*(h+hc)],[wn,h],3)
	bfacemask!(grid,[wn,0.0],[wn,0.5*(h-hc)],3)
	bfacemask!(grid,[0.0,0.0],[wn,0.0],3)
	bfacemask!(grid,[-wp,0.0],[0.0,0.0],3)

	# Interior facets
	bfacemask!(grid,[0.0,0.0],[0.0,h],4)

	# Cell regions
	cellmask!(grid,[-wp,0.0],[0.0,h],1)
	cellmask!(grid,[0.0,0.0],[wn,h],2)

	# Boundary and region types
	b_types = Dict(1=>"contact",2=>"contact",3=>"exterior",4=>"interior")
	r_types = Dict(1=>"semiconductor",2=>"semiconductor")

	# Create grid
	return (grid,b_types,r_types)

end

# 2-D planar diode with contacts at either end
function diode2d(wn,wp,h,hc;na=1e3,nd=1e3,lin=false,max_vol=0.8,
	refine_param=1e-4,yplane=nothing,bgmesh=false,x_bg=1.0,y_bg=1.0,
	max_vol_bg=2.0,dx1=0.1,dx2=0.05,dx3=0.1,ny=5,rectgrid=false)

	model = Model()

	# Material properties
	model.er = Dict(1=>11.7,2=>11.7)
	model.ni = Dict(1=>1.1e-2,2=>1.1e-2)
	model.ni_boundary = Dict(1=>1.1e-2,2=>1.1e-2,3=>1.1e-2,4=>1.1e-2)
	model.tn = Dict(1=>1e-10,2=>1e-10)
	model.tp = Dict(1=>1e-10,2=>1e-10)

	# Optional cut plane
	if yplane!=nothing
		merge!(model.ni_boundary,Dict(5=>1.1e-2))
	end

	if rectgrid
		grid,b_types,r_types = diode2d_grid_rect(
			wn,wp,h,hc,dx1,dx2,dx3,ny)
	else

		# Optional background mesh
		if bgmesh
			b,b_types,r_types = diode2d_grid(
				wn,wp,h,hc,max_vol,refine_param,yplane,
				builder=true)
			grid,b_types,r_types = bgmesh!(b,b_types,r_types,6,3,
				-wp-x_bg,wn+x_bg,-y_bg,h+y_bg,max_vol_bg)
			merge!(model.er,Dict(3=>1.0))
			merge!(model.ni_boundary,Dict(6=>1.1e-2))
		else
			grid,b_types,r_types = diode2d_grid(
				wn,wp,h,hc,max_vol,refine_param,yplane)
		end

	end

	# Linear or abrupt doping profile
	function doping_lin(x,y,reg;boundary=false)
		return (wp*nd-wn*na+(nd+na)*x)/(wp+wn)
	end
	function doping_abrupt(x,y,reg;boundary=false)
		a = 1.0e2
		b = (nd-na)/(nd+na)
		return (1/2*(nd+na)*(tanh(a*x-atanh(b))+b))
	end
	if lin
		model.doping = doping_lin
	else
		model.doping = doping_abrupt
	end

	return (Device(model,grid,b_types,r_types))

end


################################################################################
# Thyristors
################################################################################

# Thyristor from Tso 1991, "Pseudo arc-length continuation method for multiple
# solutions in one-dimensional steady-state semiconductor device simulation"
function thy1d_grid(x1,x2,x3,x4,x5,h1,h2,h3,h4,h5)

	# Create grid refined about junction center
	grid_1 = geomspace(x1,x2,h1,h2)
	grid_2 = geomspace(x2,x3,h2,h3)[2: end]
	grid_3 = geomspace(x3,x4,h3,h4)[2: end]
	grid_4 = geomspace(x4,x5,h4,h5)[2: end]
	grid_list = vcat(grid_1,grid_2,grid_3,grid_4)
	grid = simplexgrid(grid_list)

	# Create P and N-type regions
	cellmask!(grid,[x1],[x2],1)
	cellmask!(grid,[x2],[x3],2)
	cellmask!(grid,[x3],[x4],3)
	cellmask!(grid,[x4],[x5],4)
	bfacemask!(grid,[x2],[x2],3)
	bfacemask!(grid,[x3],[x3],3)
	bfacemask!(grid,[x4],[x4],3)

	# Boundary and region types
	b_types = Dict(
		1=>"contact",
		2=>"contact",
		3=>"interior"
	)
	r_types = Dict(
		1=>"semiconductor",
		2=>"semiconductor",
		3=>"semiconductor",
		4=>"semiconductor"
	)

	return (grid,b_types,r_types)

end

# Thyristor from Tso 1991, "Pseudo arc-length continuation method for multiple
# solutions in one-dimensional steady-state semiconductor device simulation"
function thy1d(;
	x1=-25.0,
	x2=-10.0,
	x3=15.0,
	x4=22.5,
	x5=25.0,
	h1=0.25,
	h2=0.25,
	h3=0.25,
	h4=0.25,
	h5=0.25,
	n1=1.0e5,
	n2=1.0e2,
	n3=1.0e3,
	n4=1.0e6
)

	grid,b_types,r_types = thy1d_grid(x1,x2,x3,x4,x5,h1,h2,h3,h4,h5)
	model = Model()

	# Material properties
	er_si = 11.7
	ni_si = 1.1e-2
	tn_si = 1.0e-10
	tp_si = 1.0e-10
	model.er = Dict(1=>er_si,2=>er_si,3=>er_si,4=>er_si)
	model.ni = Dict(1=>ni_si,2=>ni_si,3=>ni_si,4=>ni_si)
	model.ni_boundary = Dict(1=>ni_si,2=>ni_si,3=>ni_si)
	model.tn = Dict(1=>tn_si,2=>tn_si,3=>tn_si,4=>tn_si)
	model.tp = Dict(1=>tp_si,2=>tp_si,3=>tp_si,4=>tp_si)

	# Abrupt doping profile, no doping between junctions
	function doping_abrupt(x,y,reg;boundary=false)
		if !boundary
			if reg==1
				return (-n1)
			elseif reg==2
				return (n2)
			elseif reg==3
				return (-n3)
			elseif reg==4
				return (n4)
			else
				return (0.0)
			end
		else
			if reg==1
				return (-n1)
			elseif reg==2
				return (n4)
			else
				return (0.0)
			end
		end
	end
	model.doping = doping_abrupt

	d = Device(model,grid,b_types,r_types)
	return (d)

end


################################################################################
# BJTs
################################################################################

# Planar NPN with top base contact
function npn1_grid(we,wb,wc,h,hce,hcc,wcb,max_vol,refine_param,yplane;
	builder=false)

	b = SimplexGridBuilder(Generator=Triangulate)

	# Device dimensions
	w = we+wb+wc

	# Grid points
	p1 = point!(b,0,0)
	p2 = point!(b,0,(h-hce)/2)
	p3 = point!(b,0,(h+hce)/2)
	p4 = point!(b,0,h)
	p5 = point!(b,we,h)
	p6 = point!(b,we+(wb-wcb)/2,h)
	p7 = point!(b,we+(wb+wcb)/2,h)
	p8 = point!(b,we+wb,h)
	p9 = point!(b,w,h)
	p10 = point!(b,w,(h+hcc)/2)
	p11 = point!(b,w,(h-hcc)/2)
	p12 = point!(b,w,0)
	p13 = point!(b,we+wb,0)
	p14 = point!(b,we,0)

	# Emitter contact
	facetregion!(b,1)
	facet!(b,p2,p3)

	# Base contact
	facetregion!(b,2)
	facet!(b,p6,p7)

	# Collector contact
	facetregion!(b,3)
	facet!(b,p10,p11)

	# Non-contact facets
	facetregion!(b,4)
	facet!(b,p1,p2)
	facet!(b,p3,p4)
	facet!(b,p4,p5)
	facet!(b,p5,p6)
	facet!(b,p7,p8)
	facet!(b,p8,p9)
	facet!(b,p9,p10)
	facet!(b,p11,p12)
	facet!(b,p12,p13)
	facet!(b,p13,p14)
	facet!(b,p14,p1)

	# Interior facets
	facetregion!(b,5)
	facet!(b,p5,p14)
	facet!(b,p8,p13)

	# Boundary and region types
	b_types = Dict(1=>"contact",2=>"contact",3=>"contact",4=>"exterior",
		5=>"interior")
	r_types = Dict(1=>"semiconductor",2=>"semiconductor",3=>"semiconductor")

	# Optional horizontal cut plane
	if yplane!=nothing

		# Create facet for the cut plane
		p15 = point!(b,0,yplane)
		p16 = point!(b,w,yplane)
		facetregion!(b,6)
		facet!(b,p15,p16)
		merge!(b_types,Dict(6=>"cutplane"))

		# Divide into E, B and C regions, including cut plane
		cellregion!(b,1)
		regionpoint!(b,we/2,yplane/2)
		regionpoint!(b,we/2,(h+yplane)/2)
		cellregion!(b,2)
		regionpoint!(b,we+wb/2,yplane/2)
		regionpoint!(b,we+wb/2,(h+yplane)/2)
		cellregion!(b,3)
		regionpoint!(b,we+wb+wc/2,yplane/2)
		regionpoint!(b,we+wb+wc/2,(h+yplane)/2)

	else

		# Divide into E, B and C regions without cut plane
		cellregion!(b,1)
		regionpoint!(b,we/2,h/2)
		cellregion!(b,2)
		regionpoint!(b,we+wb/2,h/2)
		cellregion!(b,3)
		regionpoint!(b,we+wb+wc/2,h/2)

	end

	# Refine about both junctions
	function unsuitable(x1,y1,x2,y2,x3,y3,area)
		xcenter = (x1+x2+x3)/3
		dist = min(abs(xcenter-we),abs(xcenter-(we+wb)))
		if area > refine_param*dist
			return 1
		else
			return 0
		end
	end

	# Create grid
	options!(b,maxvolume=max_vol)
	options!(b,unsuitable=unsuitable)
	if !builder
		return (simplexgrid(b),b_types,r_types)
	else
		return (b,b_types,r_types)
	end

end

# Planar NPN with top base contact, rectangular grid
function npn1_grid_rect(we,wb,wc,h,hce,hcc,wcb,dx1,dx2,dx3,dx4,ny)

	# Create grid arrays
	x1 = geomspace(0.0,we,dx1,dx2)
	x2 = geomspace(we,we+wb,dx2,dx3)
	x3 = geomspace(we+wb,we+wb+wc,dx3,dx4)
	x = vcat(x1,x2[2: end-1],x3)
	y = range(0.0,h,length=ny)
	grid = simplexgrid(x,y)

	# Contact regions
	bfacemask!(grid,[0.0,0.5*(h-hce)],[0.0,0.5*(h+hce)],1)
	bfacemask!(grid,[we+0.5*(wb-wcb),h],[we+0.5*(wb+wcb),h],2)
	bfacemask!(grid,[we+wb+wc,0.5*(h-hcc)],[we+wb+wc,0.5*(h+hcc)],3)

	# Exterior non-contact facets
	bfacemask!(grid,[0.0,0.0],[0.0,0.5*(h-hce)],4)
	bfacemask!(grid,[0.0,0.5*(h+hce)],[0.0,h],4)
	bfacemask!(grid,[0.0,h],[we+0.5*(wb-wcb),h],4)
	bfacemask!(grid,[we+0.5*(wb+wcb),h],[we+wb+wc,h],4)
	bfacemask!(grid,[we+wb+wc,0.5*(h+hcc)],[we+wb+wc,h],4)
	bfacemask!(grid,[we+wb+wc,0.0],[we+wb+wc,0.5*(h-hcc)],4)
	bfacemask!(grid,[0.0,0.0],[we+wb+wc,0.0],4)

	# Interior facets
	bfacemask!(grid,[we,0.0],[we,h],5)
	bfacemask!(grid,[we+wb,0.0],[we+wb,h],5)

	# Cell regions
	cellmask!(grid,[0.0,0.0],[we,h],1)
	cellmask!(grid,[we,0.0],[we+wb,h],2)
	cellmask!(grid,[we+wb,0.0],[we+wb+wc,h],3)

	# Boundary and region types
	b_types = Dict(1=>"contact",2=>"contact",3=>"contact",4=>"exterior",
		5=>"interior")
	r_types = Dict(1=>"semiconductor",2=>"semiconductor",3=>"semiconductor")

	return (grid,b_types,r_types)

end

# Planar NPN with top base contact
function npn1(we,wb,wc,h,hce,hcc,wcb;ne=1e6,nb=1e4,nc=1e3,doping=nothing,
	k1=10.0,k2=10.0,max_vol=2.0,plot_doping=false,refine_param=0.3,
	yplane=nothing,dx1=1.0,dx2=0.1,dx3=0.1,dx4=2.0,ny=6,rectgrid=false)

	if rectgrid
		grid,b_types,r_types = npn1_grid_rect(
			we,wb,wc,h,hce,hcc,wcb,dx1,dx2,dx3,dx4,ny)
	else
		grid,b_types,r_types = npn1_grid(we,wb,wc,h,hce,hcc,wcb,max_vol,
			refine_param,yplane)
	end
	model = Model()

	# Material properties
	model.er = Dict(1=>11.7,2=>11.7,3=>11.7)
	model.ni = Dict(1=>1.1e-2,2=>1.1e-2,3=>1.1e-2)
	model.ni_boundary = Dict(
		1=>1.1e-2,2=>1.1e-2,3=>1.1e-2,4=>1.1e-2,5=>1.1e-2)
	model.tn = Dict(1=>1e-7,2=>1e-6,3=>1e-7)
	model.tp = Dict(1=>1e-7,2=>1e-6,3=>1e-7)

	# Optional cut plane
	if yplane!=nothing
		merge!(model.ni_boundary,Dict(6=>1.1e-2))
	end

	# Select doping profile
	if doping=="smooth"

		# Smoothed abrupt doping profile
		g1(x,y) = 1/2*(ne-(ne+nb)*tanh(k1*(x-we)
			+atanh((ne-nb)/(ne+nb))))
		g2(x,y) = 1/2*(nc+(nc+nb)*tanh(k2*(x-(we+wb))
			-atanh((nc-nb)/(nc+nb))))
		model.doping = (x,y,reg;boundary=false) -> g1(x,y)+g2(x,y)
	
	elseif doping=="expbase"

		function g_expbase(x,y,reg;boundary=false)
			if x<we
				return (ne)
			elseif we<x<we+wb
				nb_1 = nb^(1.0-(x-we)/wb)
				nb_2 = nc^((x-we)/wb)
				return (-nb_1*nb_2)
			elseif x>we+wb
				return (nc)
			else
				return (0)
			end
		end
		model.doping = g_expbase

	elseif doping=="linbase"

		function g_linbase(x,y,reg;boundary=false)
			if x<we
				return (ne)
			elseif we<x<we+wb
				a = (nc-nb)/wb
				b = 1/wb*(nb*(wb+we)-nc*we)
				return (-a*x-b)
			elseif x>we+wb
				return (nc)
			else
				return (0)
			end
		end
		model.doping = g_linbase

	else

		# Abrupt doping profile
		function g_abrupt(x,y,reg;boundary=false)
			if x<we
				return (ne)
			elseif we<x<we+wb
				return (-nb)
			elseif x>we+wb
				return (nc)
			else
				return (0)
			end
		end
		model.doping = g_abrupt

	end

	# Test plot of doping
	if plot_doping
		xtest = range(0,we+wb+wc,length=1001)
		ytest = zeros(1001)
		semilogy(xtest,abs.(model.doping.(xtest,ytest,0)),c="k")
	end

	return (Device(model,grid,b_types,r_types))

end

# Vertical NPN (similar to Padre BJT lab "discrete NPN")
function npn2_grid(x1,x2,x3,x4,y1,y2,y3,c1,c2,r,max_vol,refine_param,N_arc,
	xplane;builder=false)

	b = SimplexGridBuilder(Generator=Triangulate)

	# Device dimensions
	w = x1+x2
	h = y1+y2+y3

	# Grid points
	p1 = point!(b,0,0)
	p2 = point!(b,0,y1)
	p3 = point!(b,0,y1+y2)
	p4 = point!(b,0,h)
	p5 = point!(b,x3-c1/2,h)
	p6 = point!(b,x3+c1/2,h)
	p7 = point!(b,x1,h)
	p8 = point!(b,w-x4-c2/2,h)
	p9 = point!(b,w-x4+c2/2,h)
	p10 = point!(b,w,h)
	p11 = point!(b,w,y1)
	p12 = point!(b,w,0)
	p13 = point!(b,x1,y1+y2+r)
	p14 = point!(b,x1-r,y1+y2)

	# Emitter contact
	facetregion!(b,1)
	facet!(b,p5,p6)

	# Base contact
	facetregion!(b,2)
	facet!(b,p8,p9)

	# Collector contact
	facetregion!(b,3)
	facet!(b,p1,p12)

	# Non-contact facets
	facetregion!(b,4)
	facet!(b,p1,p2)
	facet!(b,p2,p3)
	facet!(b,p3,p4)
	facet!(b,p4,p5)
	facet!(b,p6,p7)
	facet!(b,p7,p8)
	facet!(b,p9,p10)
	facet!(b,p10,p11)
	facet!(b,p11,p12)

	# Interior facets
	facetregion!(b,5)
	facet!(b,p2,p11)
	facet!(b,p3,p14)
	facet!(b,p7,p13)
	arc!(b,r,(x1-r,y1+y2+r),-90,0,p14,p13,N=N_arc)

	# Boundary and region types
	b_types = Dict(1=>"contact",2=>"contact",3=>"contact",4=>"exterior",
		5=>"interior")
	r_types = Dict(1=>"semiconductor",2=>"semiconductor",3=>"semiconductor")

	# Optional vertical cut plane
	if xplane!=nothing

		# Create facet for the cut plane
		p15 = point!(b,xplane,0)
		p16 = point!(b,xplane,h)
		facetregion!(b,6)
		facet!(b,p15,p16)
		merge!(b_types,Dict(6=>"cutplane"))

		# Divide into E, B and C regions, including cut plane
		cellregion!(b,1)
		regionpoint!(b,xplane/2,y1+y2+y3/2)
		regionpoint!(b,(xplane+x1)/2,y1+y2+(y3+r)/2)
		cellregion!(b,2)
		regionpoint!(b,xplane/2,y1+y2/2)
		regionpoint!(b,x1+x2/2,y1+y2/2)
		cellregion!(b,3)
		regionpoint!(b,xplane/2,y1/2)
		regionpoint!(b,x1+x2/2,y1/2)

	else

		# Divide into E, B and C regions without cut plane
		cellregion!(b,1)
		regionpoint!(b,x1-r,y1+y2+r)
		cellregion!(b,2)
		regionpoint!(b,x1+x2/2,y1+y2/2)
		cellregion!(b,3)
		regionpoint!(b,(x1+x2)/2,y1/2)

	end

	# Refine about both junctions
	function unsuitable(x1u,y1u,x2u,y2u,x3u,y3u,area)

		# Compute barycenter of triangle
		xcenter = (x1u+x2u+x3u)/3
		ycenter = (y1u+y2u+y3u)/3

		# Distance from B-C junction
		dist1 = abs(ycenter-y1)

		# Distance from bottom of emitter
		if xcenter>x1-r
			dist2 = Inf
		else
			dist2 = abs(ycenter-(y1+y2))
		end

		# Distance from right edge of emitter
		if ycenter<y1+y2+r
			dist3 = Inf
		else
			dist3 = abs(xcenter-x1)
		end

		# Distance from arc in emitter
		x_arc = xcenter-(x1-r)
		y_arc = ycenter-(y1+y2+r)
		if x_arc > 0 && y_arc < 0
			dist4 = abs(sqrt(x_arc^2+y_arc^2)-r)
		else
			dist4 = Inf
		end

		# Overall refinement distance
		dist = min(dist1,dist2,dist3,dist4)
		if area > refine_param*dist
			return (1)
		else
			return (0)
		end

	end

	# Create grid
	options!(b,maxvolume=max_vol)
	options!(b,unsuitable=unsuitable)
	if !builder
		return (simplexgrid(b),b_types,r_types)
	else
		return (b,b_types,r_types)
	end

end

# Vertical NPN (similar to Padre BJT lab "discrete NPN")
function npn2(x1,x2,x3,x4,y1,y2,y3,c1,c2,r;ne=1e6,nb=1e4,nc=1e3,max_vol=2.0,
	refine_param=0.3,N_arc=10,xplane=nothing)

	grid,b_types,r_types = npn2_grid(x1,x2,x3,x4,y1,y2,y3,c1,c2,r,
		max_vol,refine_param,N_arc,xplane)
	model = Model()

	# Material properties
	model.er = Dict(1=>11.7,2=>11.7,3=>11.7)
	model.ni = Dict(1=>1.1e-2,2=>1.1e-2,3=>1.1e-2)
	model.ni_boundary = Dict(
		1=>1.1e-2,2=>1.1e-2,3=>1.1e-2,4=>1.1e-2,5=>1.1e-2)
	model.tn = Dict(1=>1e-7,2=>1e-6,3=>1e-7)
	model.tp = Dict(1=>1e-7,2=>1e-6,3=>1e-7)

	# Optional cut plane
	if xplane!=nothing
		merge!(model.ni_boundary,Dict(6=>1.1e-2))
	end

	# Abrupt doping profile
	function doping(x,y,reg;boundary=false)
		if !boundary
			if reg==1
				return (ne)
			elseif reg==2
				return (-nb)
			elseif reg==3
				return (nc)
			end
		else
			if reg==1
				return (ne)
			elseif reg==2
				return (-nb)
			elseif reg==3
				return (nc)
			else
				return (nc)
			end
		end
	end
	model.doping = doping

	return (Device(model,grid,b_types,r_types))

end

# Planar MOSFET
function mos1_grid(l,ls,h,hs,tox,c1,c2,N_arc,max_vol,builder,
	ld,hd,rs,rd,xs,xd,tch,c3,lox,max_vol_ox,max_vol_s,max_vol_d,max_vol_b,
	max_vol_ch,refine_param)

	# Default parameter values
	ld = ld==nothing ? ls : ld
	hd = hd==nothing ? hs : hd
	rs = rs==nothing ? 0.25*hs : rs
	rd = rd==nothing ? 0.25*hd : rd
	xs = xs==nothing ? 0.5*ls : xs
	xd = xd==nothing ? 0.5*ld : xd
	tch = tch==nothing ? 0.5*min(hs,hd) : tch
	c3 = c3==nothing ? c2 : c3
	lox = lox==nothing ? 1.1*c1 : lox
	max_vol_ox = max_vol_ox==nothing ? max_vol : max_vol_ox
	max_vol_s = max_vol_s==nothing ? max_vol : max_vol_s
	max_vol_d = max_vol_d==nothing ? max_vol : max_vol_d
	max_vol_b = max_vol_b==nothing ? 10*max_vol : max_vol_b
	max_vol_ch = max_vol_ch==nothing ? max_vol : max_vol_ch
	refine_param = refine_param==nothing ? 0.005 : refine_param

	b = SimplexGridBuilder(Generator=Triangulate)

	# Device dimensions
	ltot = l+ls+ld

	# Facet region numbers
	facet_g = 1
	facet_s = 2
	facet_d = 3
	facet_b = 4
	facet_ext = 5
	facet_int = 6

	# Cell region numbers
	cell_ox = 1
	cell_s = 2
	cell_d = 3
	cell_ch = 4
	cell_b = 5

	# Grid points
	p1 = point!(b,0,0)
	p2 = point!(b,0,h-hs)
	p3 = point!(b,0,h)
	p4 = point!(b,xs-c2/2,h)
	p5 = point!(b,xs+c2/2,h)
	p6 = point!(b,(ltot-lox)/2,h)
	p7 = point!(b,ls,h)
	p8 = point!(b,ltot-ld,h)
	p9 = point!(b,(ltot+lox)/2,h)
	p10 = point!(b,ltot-xd-c2/2,h)
	p11 = point!(b,ltot-xd+c2/2,h)
	p12 = point!(b,ltot,h)
	p13 = point!(b,ltot,h-hd)
	p14 = point!(b,ltot,0)
	p15 = point!(b,ls-rs,h-hs)
	p16 = point!(b,ls,h-hs+rs)
	p17 = point!(b,ls,h-tch)
	p18 = point!(b,ltot-ld+rd,h-hd)
	p19 = point!(b,ltot-ld,h-hd+rd)
	p20 = point!(b,ltot-ld,h-tch)
	p21 = point!(b,(ltot-lox)/2,h+tox)
	p22 = point!(b,(ltot-c1)/2,h+tox)
	p23 = point!(b,(ltot+c1)/2,h+tox)
	p24 = point!(b,(ltot+lox)/2,h+tox)

	# Contact facets
	facetregion!(b,facet_g)
	facet!(b,p22,p23)
	facetregion!(b,facet_s)
	facet!(b,p4,p5)
	facetregion!(b,facet_d)
	facet!(b,p10,p11)
	facetregion!(b,facet_b)
	facet!(b,p1,p14)

	# Exterior facets
	facetregion!(b,facet_ext)
	facet!(b,p1,p2)
	facet!(b,p2,p3)
	facet!(b,p3,p4)
	facet!(b,p5,p6)
	facet!(b,p6,p21)
	facet!(b,p21,p22)
	facet!(b,p23,p24)
	facet!(b,p9,p24)
	facet!(b,p9,p10)
	facet!(b,p11,p12)
	facet!(b,p12,p13)
	facet!(b,p13,p14)

	# Interior facets
	facetregion!(b,facet_int)
	facet!(b,p2,p15)
	facet!(b,p6,p7)
	facet!(b,p7,p17)
	facet!(b,p7,p8)
	facet!(b,p8,p9)
	facet!(b,p8,p20)
	facet!(b,p13,p18)
	facet!(b,p16,p17)
	facet!(b,p17,p20)
	facet!(b,p19,p20)
	arc!(b,rs,(ls-rs,h-hs+rs),-90,0,p15,p16,N=N_arc)
	arc!(b,rd,(ltot-ld+rd,h-hd+rd),-180,-90,p19,p18,N=N_arc)

	# Cell regions
	cellregion!(b,cell_ox)
	maxvolume!(b,max_vol_ox)
	regionpoint!(b,ls+l/2,h+tox/2)
	cellregion!(b,cell_s)
	maxvolume!(b,max_vol_s)
	regionpoint!(b,ls/2,h-hs/2)
	cellregion!(b,cell_d)
	maxvolume!(b,max_vol_d)
	regionpoint!(b,ltot-ld/2,h-hd/2)
	cellregion!(b,cell_b)
	maxvolume!(b,max_vol_b)
	regionpoint!(b,ls+l/2,(h-min(hs,hd))/2)
	cellregion!(b,cell_ch)
	maxvolume!(b,max_vol_ch)
	regionpoint!(b,ls+l/2,h-tch/2)

	# Boundary and region types
	b_types = Dict(
		facet_g=>"gate",
		facet_s=>"contact",
		facet_d=>"contact",
		facet_b=>"contact",
		facet_ext=>"exterior",
		facet_int=>"interior"
	)
	r_types = Dict(
		cell_ox=>"insulator",
		cell_s=>"semiconductor",
		cell_d=>"semiconductor",
		cell_b=>"semiconductor",
		cell_ch=>"semiconductor"
	)

	# Refine about the drain and source diffusion
	function unsuitable(x1,y1,x2,y2,x3,y3,area)

		# Compute barycenter of triangle
		xcenter = (x1+x2+x3)/3
		ycenter = (y1+y2+y3)/3

		# Distance from bottom of source
		if 0<xcenter<ls-rs
			dist1 = abs(ycenter-(h-hs))
		else
			dist1 = Inf
		end

		# Distance from right side of source
		if h-hs+rs<ycenter<h
			dist2 = abs(xcenter-ls)
		else
			dist2 = Inf
		end

		# Distance from arc in source
		x_arc_s = xcenter-(ls-rs)
		y_arc_s = ycenter-(h-hs+rs)
		if x_arc_s>0 && y_arc_s<0
			dist3 = abs(sqrt(x_arc_s^2+y_arc_s^2)-rs)
		else
			dist3 = Inf
		end

		# Distance from bottom of drain
		if ltot-ld+rd<xcenter<ltot
			dist4 = abs(ycenter-(h-hd))
		else
			dist4 = Inf
		end

		# Distance from left side of drain
		if h-hd+rd<ycenter<h
			dist5 = abs(xcenter-(ltot-ld))
		else
			dist5 = Inf
		end

		# Distance from arc in drain
		x_arc_d = xcenter-(ltot-ld+rd)
		y_arc_d = ycenter-(h-hd+rd)
		if x_arc_d<0 && y_arc_d<0
			dist6 = abs(sqrt(x_arc_d^2+y_arc_d^2)-rd)
		else
			dist6 = Inf
		end

		# Distance from channel
		if ls<xcenter<ltot-ld && ycenter<h
			dist7 = h-ycenter
		else
			dist7 = Inf
		end

		# Overall refinement distance
		dist = min(dist1,dist2,dist3,dist4,dist5,dist6,dist7)
		return (area>refine_param*dist ? 1 : 0)
	
	# function unsuitable
	end

	# Create grid
	options!(b,unsuitable=unsuitable)
	if !builder
		return (simplexgrid(b),b_types,r_types)
	else
		return (b,b_types,r_types)
	end

# function mos1_grid
end

# Planar MOSFET
function mos1(l,ls,h,hs,tox,c1,c2;N_arc=10,max_vol=1e-4,builder=false,
	ld=nothing,
	hd=nothing,
	rs=nothing,
	rd=nothing,
	xs=nothing,
	xd=nothing,
	tch=nothing,
	c3=nothing,
	lox=nothing,
	er_si=nothing,
	er_ox=nothing,
	ni_si=nothing,
	tn_si=nothing,
	tp_si=nothing,
	ew_ox=nothing,
	n_s=nothing,
	n_d=nothing,
	n_ch=nothing,
	n_b=nothing,
	max_vol_ox=nothing,
	max_vol_s=nothing,
	max_vol_d=nothing,
	max_vol_b=nothing,
	max_vol_ch=nothing,
	refine_param=nothing
)

	# Create grid and device model
	grid,b_types,r_types = mos1_grid(l,ls,h,hs,tox,c1,c2,N_arc,max_vol,
		builder,ld,hd,rs,rd,xs,xd,tch,c3,lox,max_vol_ox,max_vol_s,
		max_vol_d,max_vol_b,max_vol_ch,refine_param)
	model = Model()

	# Material properties (can be overridden by kwargs)
	er_si = er_si==nothing ? 11.7 : er_si
	er_ox = er_ox==nothing ? 3.9 : er_ox
	ni_si = ni_si==nothing ? 1.1e-2 : ni_si
	tn_si = tn_si==nothing ? 1e-7 : tn_si
	tp_si = tp_si==nothing ? 1e-7 : tp_si

	# Value is for n+ polysilicon with ~1e20/cm3 doping (from 6.012 S98
	# notes, Lecture 6, p. 1)
	ew_ox = ew_ox==nothing ? 5.5e-1 : ew_ox

	# Doping values (can be overridden by kwargs)
	n_s = n_s==nothing ? 2e8 : n_s
	n_d = n_d==nothing ? 2e8 : n_d
	n_ch = n_ch==nothing ? -1e6 : n_ch
	n_b = n_b==nothing ? -5e4 : n_b

	# Facet region numbers
	facet_g = 1
	facet_s = 2
	facet_d = 3
	facet_b = 4
	facet_ext = 5
	facet_int = 6

	# Cell region numbers
	cell_ox = 1
	cell_s = 2
	cell_d = 3
	cell_ch = 4
	cell_b = 5

	# Map regions to property values
	model.er = Dict(cell_ox=>er_ox,cell_s=>er_si,cell_d=>er_si,
		cell_ch=>er_si,cell_b=>er_si)
	model.ni = Dict(cell_ox=>ni_si,cell_s=>ni_si,cell_d=>ni_si,
		cell_ch=>ni_si,cell_b=>ni_si)
	model.ni_boundary = Dict(facet_s=>ni_si,facet_d=>ni_si,facet_b=>ni_si)
	model.tn = Dict(cell_ox=>tn_si,cell_s=>tn_si,cell_d=>tn_si,
		cell_ch=>tn_si,cell_b=>tn_si)
	model.tp = Dict(cell_ox=>tp_si,cell_s=>tp_si,cell_d=>tp_si,
		cell_ch=>tp_si,cell_b=>tp_si)
	model.ew = Dict(facet_g=>ew_ox)

	# Abrupt doping
	function doping(x,y,reg;boundary=false)
		doping_interior = Dict(cell_s=>n_s,cell_d=>n_d,cell_ch=>n_ch,
			cell_b=>n_b)
		doping_boundary = Dict(facet_s=>n_s,facet_d=>n_d,facet_b=>n_b)
		return (boundary ? doping_boundary[reg] : doping_interior[reg])
	end
	model.doping = doping

	# Return device object
	return (Device(model,grid,b_types,r_types))

# function mos1
end

