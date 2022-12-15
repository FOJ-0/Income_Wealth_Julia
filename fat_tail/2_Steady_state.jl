function Hugget_partialEq(he,v0)
    @unpack maxit,crit = he
    
    v = similar(v0); V=similar(v0)
    V .= v0 ; v.= v0;

    for n = 1:maxit

        V .= update_V(he,v)
        Vchange = V-v
        v .= V
        distance = maximum(abs.(Vchange)) #,1
        #println("iter $n : $distance ")
        if n == maxit
            print("iter $n : $distance ")
        end 
        
        if distance<crit
            # n_con = n
            println(" Value Function Converged, Iteration = $n ")
            break
        end
    end 
    
    V, ssb, ssf ,c, k = update_V(he,V;optimal=true)
    return V, ssb, ssf ,c, k
end

function FPE(he, ssb, ssf, c, k)
    @unpack γ, r, ρ, z1, z2, λ1, Im, amin, a, w, daab, daaf, σ2, weightb, weightf, weight0, Aswitch, dab, daf, R = he

    #RECOMPUTE TRANSITION MATRIX WITH REFLECTING BARRIER AT amax
    X = -min.(ssb,0)./daab + σ2/2. *k.^2. .*weightb
    Y = -max.(ssf,0)./daaf + min.(ssb,0)./daab + σ2/2. *k.^2. .*weight0
    Z = max.(ssf,0)./daaf + σ2/2. *k.^2. .*weightf

    A1 = spdiagm(0=>Y[:,1]) + spdiagm(-1=>X[2:Im,1]) + spdiagm(1=>Z[1:Im-1,1])
    A2 = spdiagm(0=>Y[:,2]) + spdiagm(-1=>X[2:Im,2]) + spdiagm(1=>Z[1:Im-1,2])
    A1[Im,Im] = Y[Im,1] + Z[Im,1]
    A2[Im,Im] = Y[Im,2] + Z[Im,2]    

    A = [A1 spzeros(Im,Im) ; spzeros(Im,Im) A2] .+ Aswitch
     
    #WORK WITH RESCALED DENSITY \tilde{g} BELOW
    da_tilde = 0.5*(dab + daf) 
    da_tilde[1] = 0.5*daf[1]; da_tilde[Im] = 0.5*dab[Im]  
    da_stacked = vcat(da_tilde,da_tilde)
    grid_diag = spdiagm(0=>da_stacked[:])# VER EFICIENCIA [:]
    
 
    AT  = similar(A)
    AT .= transpose(A)
    b   = zeros(2*Im)
     
    #need to fix one value, otherwise matrix is singular
    i_fix = 1
    b[i_fix] = 0.1
    row = [zeros((1,i_fix-1)) 1 zeros(1,2*Im-i_fix)]
    AT[i_fix,:] = row
    
    #Solve linear system
    g_tilde = AT\b #g_tilde
    
    #rescale \tilde{g} so that it sums to 1
    g_sum = g_tilde'*ones(2*Im)
    g_tilde = g_tilde./g_sum
    
    gg = grid_diag\g_tilde #convert from \tilde{g} to g
    
    g = [gg[1:Im] gg[Im+1:2*Im]]

    check1 = transpose(g[:,1]) *ones(Im)*da_tilde
    check2 = transpose(g[:,2]) *ones(Im)*da_tilde
    
    adot = ones(Im)*[z1 z2]*w + (R-r)*k + r.*[a a] - c

    return gg, g, adot, da_tilde
end

function Logf(he, g, da_tilde)
    @unpack Im, a, γ, σ2, ρ, r, R = he
    
    G = zeros(Im, 2)
    for i in 1:Im
        G[i,1] = sum(g[1:i,1].*da_tilde[1:i])
        G[i,2] = sum(g[1:i,2].*da_tilde[1:i])        
    end   
    
    f = zeros(Im,2)
    x = log.(max.(a, 0))
    dx = zeros(Im,1)
    
    for i =2:Im
        dx[i] = x[i] - x[i-1]
        f[i,1] = (G[i,1] - G[i-1,1])/dx[i]
        f[i,2] = (G[i,2] - G[i-1,2])/dx[i]
    end
    f[1] = 0
    
    #CALCULATE THEORETICAL POWER LAW EXPONENT
    ζ = γ*(2*σ2*(ρ - r)/(R-r)^2 -1)
    
    return x, f, ζ
    
end


function update_V(he,V; optimal = false)
     # Used in the steadt state

    @unpack γ, r, R, σ2, γ, ρ, Im,amax,amin, Δ, aa, z, zz,A,c,dVf,dVb,λ1,Aswitch, w, daab, daaf, denom, dV2b, dV2f, ϕ, weightb, weight0, weightf = he
    
    # forward difference
    dVf[1:Im-1,:] = (V[2:Im,:]-V[1:Im-1,:])./(aa[2:Im,:] - aa[1:Im-1,:])
    dVf[Im,:] = (w*z .+ r.*amax .+ (R-r)^2/(γ*σ2)*amax).^(-γ)      # will never be used
    
    # backward difference
    dVb[2:Im,:] = (V[2:Im,:]-V[1:Im-1,:])./(aa[2:Im,:] - aa[1:Im-1,:])
    dVb[1,:] = (w*z .+ r.*amin).^(-γ) # state constraint boundary condition

    
    
    #second derivative: approximation only differs at amax
    dV2b[2:Im-1,:] = (daab[2:Im-1,:].*V[3:Im,:] - (daab[2:Im-1,:] + daaf[2:Im-1,:]).*V[2:Im-1,:] + daaf[2:Im-1,:].*V[1:Im-2,:])./denom[2:Im-1,:]
    dV2f[2:Im-1,:] = (daab[2:Im-1,:].*V[3:Im,:] - (daab[2:Im-1,:] + daaf[2:Im-1,:]).*V[2:Im-1,:] + daaf[2:Im-1,:].*V[1:Im-2,:])./denom[2:Im-1,:]
    dV2b[Im,:] = -γ*dVb[Im,:]/amax
    dV2f[Im,:] = -γ*dVf[Im,:]/amax
    
    
    
    # consumption and savings with forward difference
    cf = max.(dVf,10^(-10)).^(-1/γ)
    kf = max.(- dVf./dV2f.*(R-r)/σ2,0)  
    kf = min.(kf, aa .+ ϕ) 
    ssf = w*zz + (R-r).*kf + r.*aa - cf

    # println((sum(cf), sum(kf), sum(ssf))) 
    
    # consumption and savings with backward difference
    cb = max.(dVb,10^(-10)).^(-1/γ)
    kb = max.(- dVb./dV2b.*(R-r)/σ2,0)
    kb = min.(kb, aa .+ ϕ)
    ssb = w*zz + (R-r).*kb + r.*aa - cb
 
    # consumption and derivative of value function at steady state
    k0 = (kb + kf)/2 #could do something fancier here but the following seems to work well. And more importantly, it's fast
    c0 = w*zz + (R-r).*k0 + r.*aa
    dV0 = max.(c0, 10^(-10)).^(-γ);
    
    # dV_upwind makes a choice of forward or backward differences based on
    # the sign of the drift    
    If = ssf .> 10^(-12) # positive drift --> forward difference 
    Ib = ssb .< -10^(-12).*(1 .-If) # negative drift --> backward difference
    I0 = (1 .- If .-Ib)

#     # STATE CONSTRAINT at amin: USE BOUNDARY CONDITION UNLESS muf > 0:
    dV_Upwind = dVf.*If + dVb.*Ib + dV0.*I0 #important to include third term
    c[:] = max.(dV_Upwind,10^(-10)).^(-1/γ)
    u = c.^(1-γ)/(1-γ)
    k = max.(-dV_Upwind./dV2b.*(R-r)/σ2,0)
    k = min.(k, aa .+ ϕ);    
    
    X = -Ib.*ssb./daab + σ2/2. *k.^2. .*weightb
    Y = - If.*ssf./daaf + Ib.*ssb./daab + σ2/2. *k.^2. .*weight0
    Z = If.*ssf./daaf + σ2/2. *k.^2. .*weightf
  
    ξ = -amax*(R-r)^2/(2*γ*σ2)
    X[Im,:] = -min.(ssb[Im,:],0)./daab[Im,:] - ξ./daab[Im,:]
    Y[Im,:] = -max.(ssf[Im,:],0)./daaf[Im,:] + min.(ssb[Im,:],0)./daab[Im,:] + ξ./daab[Im,:]
    Z[Im,:] = max.(ssf[Im,:],0)./daaf[Im,:]
    
    
    A1 = spdiagm(0=>Y[:,1]) + spdiagm(-1=>X[2:Im,1]) + spdiagm(1=>Z[1:Im-1,1])
    A2 = spdiagm(0=>Y[:,2]) + spdiagm(-1=>X[2:Im,2]) + spdiagm(1=>Z[1:Im-1,2])
    A1[Im,Im] = Y[Im,1] + Z[Im,1]
    A2[Im,Im] = Y[Im,2] + Z[Im,2]    
    
    A = [A1 spzeros(Im,Im) ; spzeros(Im,Im) A2] .+ Aswitch
    
    B = (ρ + 1/Δ) .* sparse(1.0I, 2*Im, 2*Im) .- A

    u_stacked = vcat(u[:,1],u[:,2])
    V_stacked = vcat(V[:,1],V[:,2])

    b = u_stacked .+ V_stacked/Δ
    V_stacked = B\b    #SOLVE SYSTEM OF EQUATIONS
    V = hcat(V_stacked[1:Im],V_stacked[Im+1:2*Im]) 
    
    if optimal == true
        return V, ssb, ssf ,c, k
    else
        return V
    end
end
