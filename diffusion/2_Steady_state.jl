function Hugget_steadystate(he,v0)
    @unpack critK, maxitK, J, Im = he
    @unpack K, r, w = he #Initial values

    V = similar(v0)
    c=0.; g=0.; gg=0.; S=0.; 
    
    for i in 1:maxitK
        he    = HuggetEconomy(r=r, w=w, K=K, J=J, Im=Im)
        V, A, c = Hugget_partialEq(he,v0)
        S, gg, g = FPE(he, A, c)
        K, r, w  = update_prices(he, K, S)
        v0 .= V
        
        if abs(K-S) < critK
            println("Convergence ok, iter:", i)
            break
        end

        if i == maxitK
            println("Not converge, i = maxitK:", i)
        end
    end
    
    he = HuggetEconomy(r=r, w=w, K=K, J=J, Im=Im)
    V, A, c = Hugget_partialEq(he, V)
    S, gg, g = FPE(he, A, c)
    ss = f_ss(he, c)
    return V, c, g, gg, S, ss

end


function Hugget_partialEq(he,v0)
    @unpack Im,maxit,crit = he
    
    v = similar(v0); V=similar(v0)
    V .= v0 ; v.= v0;

    for n = 1:maxit

        V .= update_V(he,v)
        Vchange = V-v
        v .= V
        distance = maximum(abs.(Vchange)) 

        if n == maxit
            print("iter $n : $distance ")
        end 
        
        if distance<crit
            n_con = n
            break
        end
    end 
    
    V,A,c = update_V(he,V;optimal=true)
    return V,A,c
end

function FPE(he, A, c)
    @unpack Im, a, da, dz, J = he
    
    AT  = similar(A)
    AT .= transpose(A)
    b   = zeros(J*Im)

    #need to fix one value, otherwise matrix is singular
    i_fix = 1
    b[i_fix] = 0.1
    
    for j=1:Im*J
        AT[i_fix,j]=0
    end
    
    AT[i_fix,i_fix]=1    
    #Solve linear system
    gg = AT\b
    g_sum = gg'*ones(J*Im)*da*dz
    gg = gg./g_sum
    g = reshape(gg,Im,J)
    S = sum(g.*a*da*dz)

    return S, gg, g
end

function f_ss(he, c)
    @unpack w, zz, r, aa = he
    return w*zz + r.*aa - c
end


function update_prices(he, K, S)
    @unpack relax, α, δ = he
    K = relax*K +(1-relax)*S           #relaxation algorithm (to ensure convergence)
    r = α * K^(α-1) - δ #interest rates
    w = (1-α) * K^(α)     
    
    return K, r, w
end

function update_V(he,V; optimal = false)
     # Used in the steadt state

    @unpack σ, r, ρ, Im,amax,amin, Δ, da, aa,z, zz,A,c,dVaf,dVab,Aswitch, w, J, centdiag, lowdiag, updiag = he
    # forward difference
    dVaf[1:Im-1,:] = (V[2:Im,:]-V[1:Im-1,:])/da
    dVaf[Im,:] = (w*z .+ r.*amax).^(-σ)      # will never be used

    # backward difference
    dVab[2:Im,:] = (V[2:Im,:]-V[1:Im-1,:])/da;
    dVab[1,:] = (w*z .+ r.*amin).^(-σ); # state constraint boundary condition

    # consumption and savings with forward difference
    cf = dVaf.^(-1/σ);
    ssf = w*zz + r.*aa - cf

    # consumption and savings with backward difference
    cb = dVab.^(-1/σ);
    ssb = w*zz + r.*aa - cb;

    # consumption and derivative of value function at steady state
    c0 = w*zz + r.*aa;
    dV0 = c0.^(-σ)

    # dV_upwind makes a choice of forward or backward differences based on
    # the sign of the drift    
    If = ssf .> 0; # positive drift --> forward difference
    Ib = ssb .< 0; # negative drift --> backward difference
    I0 = (1 .- If .-Ib)

    # STATE CONSTRAINT at amin: USE BOUNDARY CONDITION UNLESS muf > 0:
    dV_Upwind = dVaf.*If + dVab.*Ib + dV0.*I0; #important to include third term
    
    c[:] = dV_Upwind.^(-1/σ)
    u = c.^(1-σ)/(1-σ)

    X = - min.(ssb,0)/da
    Y = - max.(ssf,0)/da + min.(ssb,0)/da
    Z = max.(ssf,0)/da

	for j=1:J
        for i=1:Im
            @inbounds centdiag[Im*(j-1)+i] = Y[i,j];
            if i<Im
                @inbounds updiag[Im*(j-1)+i+1] = Z[i,j];
                @inbounds lowdiag[Im*(j-1)+i] = X[i+1,j];
            end
        end
	end    

    
    A = spdiagm(Im*J, Im*J, 0 => centdiag) + spdiagm(Im*J, Im*J, -1 => lowdiag[1:Im*J-1]) + spdiagm(Im*J, Im*J, 1 => updiag[2:Im*J]) + Aswitch
    B = (ρ + 1/Δ) .* sparse(1.0I, J*Im, J*Im) .- A
    
    u_stacked = reshape(u,Im*J,1) #VER EFICIENCIA
    V_stacked = reshape(V,Im*J,1) #VER EFICIENCIA

    b = u_stacked .+ V_stacked/Δ
    V_stacked = B\b    #SOLVE SYSTEM OF EQUATIONS
    V = reshape(V_stacked,Im,J)

    if optimal == true
        return V,A,c
    else
        return V
    end

end
