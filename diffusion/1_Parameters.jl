function HuggetEconomy(;σ=2.0, # CRRA utility with parameter gamma
        α= 0.35,               # Production function F = K^alpha * L^(1-alpha) 
        δ= 0.1,                # Capital depreciation
        zmean = 1.0,           # mean O-U process (in levels). This parameter has to be adjusted to ensure that the mean of z (truncated gaussian) is 1.
        σ2 = (0.10)^2,         # sigma^2 O-U
        Corr = exp(-0.3),      # persistence -log(Corr)  O-U
        ρ = 0.05,              # discount rate
        K = 3.8,               # initial aggregate capital. It is important to guess a value close to the solution for the algorithm to converge
        relax = 0.99,          # relaxation parameter 
        J=40,                  # number of z points 
        zmin=0.5, zmax=1.5,    #zmin, zmax
        amin = -1,             # borrowing constraint
        amax = 30.,            # range a
        Im = 1000,             # number of a points 
        maxit=100,             #maximum number of iterations in the HJB loop
        maxitK = 100,          #maximum number of iterations in the K loop
        crit = 10^(-6.0),      #criterion HJB loop
        critK = 10^(-5.),      #criterion K loop
        Δ = 1000.,             #delta in HJB algorithm
        r = 0.,
        w = 0.,
        
        )

    the = -log(Corr);
    Var = σ2/(2*the);    
    
    #VARIABLES     
    a = collect(range(amin,stop=amax,length=Im))
    # aa = [a a]
    da = (amax-amin)/(Im-1)     
    

    z = collect(range(zmin,stop=zmax,length=J))'  
    dz = (zmax-zmin)/(J-1)
    dz2 = dz*dz  

    aa = a*ones(1,J)
    zz = ones(Im,1)*z 

    μ = the .*(zmean .- z)        #DRIFT (FROM ITO'S LEMMA)
    s2 = σ2.*ones(1,J)        #VARIANCE (FROM ITO'S LEMMA)
 
    r = α * K^(α-1) - δ
    w = (1-α) * K^(α)

    #Matrix dim
    χ =  - min.(μ,0) ./dz + s2/(2*dz2);
    yy =  min.(μ,0)./dz - max.(μ,0)./dz - s2/dz2;
    ζ = max.(μ,0)./dz + s2/(2*dz2);
    
    
    
    updiag = zeros(Im*J);
    centdiag = zeros(Im*J);
    lowdiag = zeros(Im*J);

    for i = 1:Im
        centdiag[i]=χ[1]+yy[1];
        for j = 1:J-2
            centdiag[Im*j+i] = yy[j+1];
            lowdiag[Im*(j-1)+i] = χ[j+1];
            updiag[Im*j+i] = ζ[j];
        end
        centdiag[(J-1)*Im+i] = yy[J]+ζ[J];
        updiag[(J-1)*Im+i] = ζ[J-1];
        lowdiag[Im*(J-2)+i] = χ[J];
    end    
    
    Aswitch = spdiagm(Im*J, Im*J, 0 => centdiag) + spdiagm(Im*J, Im*J, -Im => lowdiag[1:Im*J-Im]) + spdiagm(Im*J, Im*J, Im => updiag[Im+1:Im*J])
    
    #Clear updiag and lowdiag
    for j = 1:J
       updiag[Im*(j-1)+1] = 0;
       lowdiag[Im*j] = 0;
    end    

    dist =zeros(maxit)
    v0 = zeros(Im,J)    
    
 
    A = sparse(1.0*I,J*Im, J*Im)
    c = zeros(Im,J)
    dVaf = zeros(Im,J); dVab = zeros(Im,J);
    dVzf = zeros(Im,J); dVzb = zeros(Im,J);

    
    return(σ=σ, α=α, δ=δ, zmean=zmean, σ2=σ2, Corr=Corr, ρ=ρ, K=K, relax=relax, J=J, zmin=zmin,
        zmax=zmax, amin=amin, amax=amax, Im=Im, maxit=maxit, maxitK=maxitK, crit=crit, critK=critK, 
        Δ=Δ, the=the, Var=Var, a=a, da=da, z=z, dz=dz, dz2=dz2, aa=aa, zz=zz, μ=μ, s2=s2, Aswitch=Aswitch,
        dist=dist, v0=v0, A=A, c=c, dVaf=dVaf, dVab=dVab, dVzf=dVzf, dVzb=dVzb, r=r, w=w, 
        centdiag=centdiag, lowdiag=lowdiag, updiag=updiag) 
    
end

function initial_V(he)
    @unpack σ, r, ρ,aa, zz,Im, w = he
    v0 = (w*zz + r.*aa).^(1-σ)/(1-σ)/ρ
    return v0
end