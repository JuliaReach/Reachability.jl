#=Copyright (c) 2018-2019: David Sanders and Luis Benet.
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.=#

function validated_integ(f!, qq0::AbstractArray{T,1}, δq0::IntervalBox{N,T},
        t0::T, tmax::T, orderQ::Int, orderT::Int, abstol::T;
        maxsteps::Int=500, parse_eqs::Bool=true,
        check_property::Function=(t, x)->true) where {N, T<:Real}
    println("hello mom!")
    # Set proper parameters for jet transport
    @assert N == get_numvars()
    dof = N
    # if get_order() != orderQ
    #     set_variables("δ", numvars=dof, order=orderQ)
    # end

    # Some variables
    R   = Interval{T}
    q0 = IntervalBox(qq0)
    t   = t0 + Taylor1(orderT)
    tI  = t0 + Taylor1(orderT+1)
    δq_norm = IntervalBox(Interval{T}(-1, 1), Val(N))
    q0box = q0 .+ δq_norm

    # Allocation of vectors
    # Output
    tv    = Array{T}(undef, maxsteps+1)
    xv    = Array{IntervalBox{N,T}}(undef, maxsteps+1)
    # xTMNv = Array{TaylorModelN{N,T,T}}(undef, dof, maxsteps+1)
    # Internals: jet transport integration
    x     = Array{Taylor1{TaylorN{T}}}(undef, dof)
    dx    = Array{Taylor1{TaylorN{T}}}(undef, dof)
    xaux  = Array{Taylor1{TaylorN{T}}}(undef, dof)
    x0    = Array{TaylorN{T}}(undef, dof)
    xTMN  = Array{TaylorModelN{N,T,T}}(undef, dof)
    # Internals: Taylor1{Interval{T}} integration
    xI    = Array{Taylor1{Interval{T}}}(undef, dof)
    dxI   = Array{Taylor1{Interval{T}}}(undef, dof)
    xauxI = Array{Taylor1{Interval{T}}}(undef, dof)
    x0I   = Array{Interval{T}}(undef, dof)

    # Set initial conditions
    zI = zero(R)
    Δ = zero.(q0)
    rem = Array{Interval{T}}(undef, dof)
    @inbounds for i in eachindex(x)
        qaux = normalize_taylor(qq0[i] + TaylorN(i, order=orderQ), δq0, true)
        x[i] = Taylor1( qaux, orderT)
        dx[i] = x[i]
        x0[i] = copy(qaux)
        xTMN[i] = TaylorModelN(x[i][0], zI, q0, q0box)
        #
        xI[i] = Taylor1( q0box[i], orderT+1 )
        dxI[i] = xI[i]
        x0I[i] = qaux(δq_norm)
        rem[i] = zI
    end

    # Output vectors
    @inbounds tv[1] = t0
    @inbounds xv[1] = IntervalBox( evaluate(xTMN, δq_norm) )
    # @inbounds xTMNv[:, 1] .= xTMN[:]

    # Determine if specialized jetcoeffs! method exists (built by @taylorize)
    parse_eqs = parse_eqs && (length(methods(TaylorIntegration.jetcoeffs!)) > 2)
    if parse_eqs
        try
            TaylorIntegration.jetcoeffs!(Val(f!), t, x, dx)
        catch
            parse_eqs = false
        end
    end

    # Integration
    nsteps = 1
    while t0 < tmax

        # Validated step of the integration
        δt = validated_step!(f!, t, x, dx, xaux, tI, xI, dxI, xauxI,
            t0, tmax, x0, x0I, xTMN, xv, rem, δq_norm,
            q0, q0box, nsteps, orderT, abstol, parse_eqs, check_property)

        # New initial conditions and time
        nsteps += 1
        t0 += δt
        @inbounds t[0] = t0
        @inbounds tI[0] = t0
        @inbounds tv[nsteps] = t0
        @inbounds for i in eachindex(x)
            aux = x[i](δt)
            x[i]  = Taylor1( aux, orderT )
            dx[i] = Taylor1( zero(aux), orderT )
        end
        # @show(IntervalBox(rem))
        # @inbounds xTMNv[:, nsteps] .= xTMN[:]

        # println(nsteps, "\t", t0, "\t", remainder.(xTMN[:]), "\t", diam(Δ))
        if nsteps > maxsteps
            @info("""
            Maximum number of integration steps reached; exiting.
            """)
            break
        end

    end

    return view(tv,1:nsteps), view(xv,1:nsteps)#, view(xTMNv, 1:nsteps, :)
end

