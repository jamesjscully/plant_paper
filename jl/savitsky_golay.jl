struct SavitzkyGolayFilter{M,N} end
@generated function (::SavitzkyGolayFilter{M,N})(data::AbstractVector{T}) where {M, N, T}
            #Create Jacobian matrix
            J = zeros(2M+1, N+1)
            for i=1:2M+1, j=1:N+1
                J[i, j] = (i-M-1)^(j-1)
            end
            e₁ = zeros(N+1)
            e₁[1] = 1.0

            #Compute filter coefficients
            C = J' \ e₁

            #Evaluate filter on data matrix

            To = typeof(C[1] * one(T)) #Calculate type of output
            expr = quote
                n = size(data, 1)
                smoothed = zeros($To, n)
                @inbounds for i in eachindex(smoothed)
                    smoothed[i] += $(C[M+1])*data[i]
                end
                smoothed
            end

            for j=1:M
                insert!(expr.args[6].args[3].args[2].args, 1,
                    :(if i - $j ≥ 1
                        smoothed[i] += $(C[M+1-j])*data[i-$j]
                      end)
                )
                push!(expr.args[6].args[3].args[2].args,
                    :(if i + $j ≤ n
                        smoothed[i] += $(C[M+1+j])*data[i+$j]
                      end)
                )
            end

            return expr
end
