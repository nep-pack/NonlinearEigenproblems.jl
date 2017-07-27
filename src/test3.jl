function f(A::Array{Float64,2}, b::Vector{Float64})
          x = zeros(3);
          m = 0.0;
          al,ac = size(A)
          for k= 1:(al-1)
              #println("valor de k ",k)
              for i = (k+1):al
                  #println("valor de i ",i)
                  m = A[i,k]/(A[k,k])
                  A[i,k] = 0
                  for j=(k+1):al
                      #println("valor de j",j)
                      A[i,j] = A[i,j] - m*A[k,j]
                      b[i]= b[i] - m*b[k]
                  end
              end
          end
          x[al] = b[al]/(A[al,al])
          for k = (al-1):-1:1
              begin
              s = 0;
              for j = (k+1):al
                  s = s+A[k,j]*x[j]
              end
              x[k]=(b[k]-s)/A[k,k]
          end
      end
      return x
end
