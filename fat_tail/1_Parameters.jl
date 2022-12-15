function Economy(;γ=2.0, 
        ρ = 0.05, 
        z1=0.01, z2=3*z1,         
        λ1=0.5, λ2=0.5, 
        w= 3.0,
        R = 0.051,
        r0=0.041, 
        ϕ = 0.3,
        ζ = 1.5,
        Im = 5000,
        amin = -ϕ, amax = 1000., 
        coeff = 5, power = 10,
        maxit=100, 
        crit = 10^(-6),
        Δ = 1000.,        
        )
    
    r = r0
    z  = [z1 z2]
    λ  = [λ1 λ2]

    σ2 = (ζ/γ + 1)/2*(R-r)^2/(ρ - r) #pick sig2 so that zeta=1.5
    σ = σ2^(1/2)


    x = collect(range(0,stop=1,length=Im))
    xx = x + coeff*x.^power
    xmax = maximum(xx); xmin = minimum(xx)
    a = (amax-amin)/(xmax - xmin) .*xx .+ amin
    daf = ones(Im,1); dab = ones(Im,1)
    daf[1:Im-1] = a[2:Im] - a[1:Im-1]
    dab[2:Im]   = a[2:Im] - a[1:Im-1]
    daf[Im]     = daf[Im-1]
    dab[1]      = dab[2]
    
    aa = [a a]
    daaf = daf*ones(1,2)
    daab = dab*ones(1,2) 
    
    #objects for approximation of second derivatives
    denom = 0.5*(daaf + daab).*daab.*daaf
    weightf = daab./denom
    weight0 = -(daab + daaf)./denom
    weightb = daaf./denom
    
    zz = ones(Im)*z
    kk = ones(Im,2)
        
    dist =zeros(maxit)
    v0 = zeros(Im,2)
    c = zeros(Im,2)
    dVf = zeros(Im,2); dVb = zeros(Im,2)      
    dV2f = zeros(Im,2); dV2b = zeros(Im,2)    
     
    Aswitch = [ -sparse(I, Im, Im)*λ[1]  sparse(I, Im, Im)*λ[1] ;
                sparse(I, Im, Im)*λ[2]   -sparse(I, Im, Im)*λ[2] ]
    
    A = sparse(1.0*I,2*Im, 2*Im)


        return(γ=γ, ρ=ρ, z1=z1, z2=z2, λ1=λ1, λ2=λ2, w=w, R=R, r0=r0, ϕ=ϕ, ζ=ζ, Im=Im, amin=amin, amax=amax,
               coeff=coeff, power=power, maxit=maxit, crit=crit, Δ=Δ, r=r, z=z, λ=λ, σ2=σ2, σ=σ, 
               a=a, daf=daf, dab=dab, aa=aa, daaf=daaf, daab=daab, denom=denom, weightf=weightf, weight0=weight0, 
               weightb=weightb, zz=zz, kk=kk, dist=dist, v0=v0, c=c, dVf=dVf, dVb=dVb, Aswitch=Aswitch, A=A,
               dV2b=dV2b, dV2f=dV2f)
end

function initial_V(he)
    @unpack γ, r, ρ,a, z,Im, w = he
    v0 = zeros(Im,2)
    v0[:,1] = ((w*z[1] .+ r.*a).^(1-γ)/(1-γ)/ρ)
    v0[:,2] = ((w*z[2] .+ r.*a).^(1-γ)/(1-γ)/ρ);
    return v0
end