
#COSES GAUSS
function dot(x,y)
    result  = 0
    for i in 1:size(x)[1]
        result += x[i] * y[i]
    end
    
    return result
end

function norm(v)
    return sqrt(dot(v,v))
end


function Gauss(A, b)
    n, m = size(A)
    Ab = hcat(A, b)
    for i in 1:n-1
        # Partial pivoting
        max_index = i
        for k in i+1:n
            if abs(Ab[k, i]) > abs(Ab[max_index, i])
                max_index = k
            end
        end
        # Swap 
        Ab[[i, max_index], :] = Ab[[max_index, i], :]
        
        # Check for singular matrix
        if Ab[i, i] == 0
            prinln("Singular matrix")
        end
        
        # Elimination step
        for j in i+1:n
            factor = Ab[j, i] / Ab[i, i]
            Ab[j, i:end] -= factor * Ab[i, i:end]
        end
    end
    
    # Back substitution
    x = zeros(n)
    x[n] = Ab[n, end] / Ab[n, n]
    for i in n-1:-1:1
        x[i] = (Ab[i, end] - dot(Ab[i, i+1:end-1], x[i+1:end])) / Ab[i, i]
    end
    
    return x
end

function identity(n)
     I = zeros(n,n)
     for i in 1:n
        I[i,i] = 1
     end
     return I
end

#POWER method
function power_method(A, max_iter, tol)
    n = size(A)[1]
    # Initial guess for the eigenvector
    v = zeros(n)
    v[1] = 1.0
    
    a_0 = 0.0
    a = 0.0
    
    for iter in 1:max_iter
        # Multiply A by the current estimate of the eigenvector
        Av = A * v
        
        # Compute the new estimate of the eigenvalue
        a = dot(v, Av)
        
        # Normalize the eigenvector
        v = Av / norm(Av)
        
        # Check for convergence
        if abs(a-a_0) < tol
            return a
        end
        
        a_0 = a
    end
    return a
end

#COSES POLINOMIS
function reformat(x,n)
    m = min(n,size(x)[1])
    y = [0.0 for i in 1:n]
    for i in 1:m
        y[i] = x[i]
    end
    return y 
end
function truncate(v,n)
    n = n+1
    v = reformat(v,n)
    for i in 2:2:n
        v[i] = 0
    end
    return [ v[i] for i in 1:n ]
end

function vec_to_pol(x)
    if(x[1]!= 0) 
        print(x[1],"+")
    end
    for i in 2: size(x)[1]
        if(x[i]!= 0) 
            print(x[i],"X^(",(i-1),")+")
        end
    end
end

function sum(x,y)
    n = size(x)[1]
    m = size(y)[1]
    k = max(n,m)
    z = [0 for i in 1:k]
    x = reformat(x,k)
    y = reformat(y,k)
    return x+y
end

function product(x,y)
    n = size(x)[1]

    m = size(y)[1]
    z = [ 0.0 for i in 1:(n+m)]
    for i in 1: n
        for j in 1:m
        z[i+j-1] = z[i+j-1] + x[i]*y[j] 
        end
    end
    return z
end

function composition(x,y)
    n = size(x)[1]
    m = size(y)[1]
    z = [x[n]]
    for i in n:-1:2
        z  = product(z,y)
        z[1] = z[1] + x[i-1]
    end
    
    return reformat(z,(n-1)*(m-1)+1)
end
function evaluate(v,a)
    return composition(v,[a])[1]
end
function derivate(x)
    n = size(x)[1]
    return [ x[i+1]*i for i in 1:n-1]
end

function power(x,n)
    z = x
    for i in  1:n-1
        z = product(z,x)
    end
    return z
end

function T(ϕ)
    a = evaluate(ϕ,1)
    x = [0.0,1.0]
    result = a * x 
    result = composition(ϕ,result)
    result = composition(ϕ,result)
    result = result * 1.0/a
    return  result
end


function DT_ai(ϕ,i,n)
    a = evaluate(ϕ,1.0)
    x = [0.0,1.0]
    p = -1*composition(ϕ,composition(ϕ,a*x))/a^2
    q = sum(composition(power(x,i),a*x), product(composition(derivate(ϕ),a*x),x))
    q = product(q,composition(derivate(ϕ),composition(ϕ,a*x)))
    q = sum(q,composition(power(x,i),composition(ϕ,a*x)))
    p = sum(p,q* 1.0/a)
    p = truncate(p,n)
    return p[2:end]
end

function DT(ϕ,n)
    DT = zeros(n,n)
    #coomputem DT començant desde i= 2 que equival a x^2 i saltan de 2 en  2
    for i in 2:2:n
        DT[:,i] = DT_ai(ϕ,i,n)
    end
    #retornem nomes la part parella
    return  DT[2:2:end,2:2:end]
end

function Newton(ϕ,n,tol)
    i = 10^10
    x = [0.0, 1.0]
    while (true) 
        DF =  - DT(ϕ,n)
        k = size(DF)[1]
        I = identity(k)
        DF = DF + I
        F = reformat(sum(ϕ,- T(ϕ)),2*k+3)
        A = -DF
        b = F[3:2:n+1]
        Δa = Gauss(A,b)
        ϕ = sum(ϕ,composition(product(Δa,x),[0.0,0.0,1.0]))
        if (norm(Δa)<tol)
            #println("Proces de Newton acabat exitosament",Δa) # serveix per veure si tot funciona
            return ϕ
        end
        i = i-1
        if (i==0)
            println("Error")
            return 0
        end
    end
end


function main()
    ϕ = [1.0,0,-1.5]
    i = 40
    open("Polinomis2.txt","w") do io
        for i in 4:2:40
            println("Estat ", i/40*100)
            #Newton per resoldre 
            ϕ = Newton(ϕ,i,10^(-15))
            #Powm es el powet method del paquet iterative solvers, retorna el vap i el vep, ens quedem nomes amb el vap [1] i agafem la part real
            μ = power_method(DT(ϕ,i),10^6,10^(-15))
            #println(i,"&",round(μ,digits=15),"\\\\") (imprimeix per posar a LaTex)
            println(io,"Polinomi de grau ",i, " ϕ = ",ϕ, " amb error: ",norm(sum(T(ϕ),-ϕ)), " constant de Feigenbaum ≈ ", μ )
        end
    end
end

main()


