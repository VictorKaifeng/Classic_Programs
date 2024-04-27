using Plots
# Elprograma esta escrit per blocs, es recomana copiar cada block a unfitxer diferent i executarlo 
# o anar executant cada bloc per separat estil Notebook

# Escriure la funció
function f(x,μ)
    return μ*x*(1-x)
end

function f_k(x,μ,k)
    y = 0
    for i in 1:k
        y = f(x,μ)
        x = y
    end
    return x
end
function partial_x(x,a,k)
    h = 10^(-8)
    return (f_k(x+h,a,k)-f_k(x,a,k))*(1.0/h)
    end
function partial_a(x,a,k)
    h = 10^(-8)
    return (f_k(x,a+h,k)-f_k(x,a,k))*(1.0/h)
    end

# Resoldre la funcio  donada a
function Newton_x(x,a,k)
    y=0
    for i in 1:10^3
        y = x - (f_k(x,a,k)-x)/(partial_x(x,a,k)-1)
        if(abs(y-x)<10^(-10))
            return y
        end
        x = y
    end
    return x
end

function solve(a,k)
        return Newton_x(0.5,a,k)
end


#Apartat 1
#fa el graph en un fitxer per gnuplot 
function graph(n1,m1,n2,m2,iterates)
    open("punts.dat","w") do io
        for i in 1:n1
            a = 1 + 3.0*i/n1
            if(i%1000==0)
                println("Estem al pas1 %",100*i/n1)
            end
            for j in 1:m1
                x = 1.0*j/m1
                if (abs(f_k(x,a,iterates)-x)<10^(-4))
                    println(io,a," ",x)
                end      
            end
        end
        for i in 1:n2
            a = 3 + 1.0*i/n2
            if(i%1000==0)
                println("Estem al pas2 amb %",100*i/n2)
            end
            for j in 1:m2
                x = 1.0*j/m2
                if (abs(f_k(x,a,iterates)-x)<10^(-4))
                    println(io,a," ",x)
                end      
            end
        end
    end
end

#graph(10^3,10^3,10^4,10^4,3^3)
#la funcio fitxer amb aquestos valors pot tardar molt per aixo es mante en comentari per evitar executarla sense voler


#Apartat 2
#trobar els intervals
function periode(x,μ,base,power)
    y = 0
    for i in 1:power
        k = base^i
        y = f_k(x,μ,k)
        if (abs(y-x)<10^(-4))
            return k
        end
    end
    return 0
end

function inici_interval(periods,k)
    k = k /2
    for i in 1:length(periods)-1
        if (periods[i][2] == periods[i+1][2])
            if (periods[i][2]>k)
                return periods[i][1]
            end
        end
    end
end

function intervals(a,b,n,base,power)
    periods = NTuple{2,Float64}[]

    for μ in range(a,b,n)
        x = solve(μ,base^power)
        k = periode(x,μ,base,power)
        push!(periods,(μ,k))
    end

    for k in 1:8
        println("El periode ",2^k," comença a: ",inici_interval(periods,2^k))
    end
    return periods
end
# Modificar valor de a i b per anar fent els diferents grafics
a = 3.6
b = 3.69
p = intervals(a,b,10^2,2,10)
p2 = [ x for x in p if x[2]!=0] # treiem els punts que no han obtingut cap periode

scatter( p2, size=(1080,720), xlims = (3.68,b), ylims =(10,100))



#APARTAT 3

#Troba les an
function secant_a(a,b,k)
    c = 0
        for i in 1:20
            x = solve(a,k)
            y = solve(b,k)
            fa = partial_x(x,a,k)+1
            fb = partial_x(y,b,k)+1
            c = (b*fa-a*fb)/(fa-fb)
            if(abs(partial_x(solve(c,k),c,k)+1)<10^(-7))
                return c
            end
            b = a
            a = c
        end
        return c
    end
function find_an(A,base)
    B = zeros(Float64,9)
    for i in 0:8
        k = base^i
        B[i+1] =  secant_a(A[i+1][1],A[i+1][2],k)
        println("Al pas ", i," tenim:", B[i+1],"\n")
    end
    return B
end

#Computa els ratios 
function ratios(A,base)
    a = find_an(A,base )
    p= zeros(9)

    println("Ratio")
    function ratio(a,p)
        for i in 3:8
            p[i-2] = (a[i-1]-a[i-2])/(a[i]-a[i-1])
            println((a[i-1]-a[i-2])/(a[i]-a[i-1]))
        end
    end

    ratio(a,p)

    println("Ratio amb Aitken")
    function Aitken(p)
        for i in 1:4
            println(p[i+2]-(p[i+2]-p[i+1])^2/(p[i+2]-2*p[i+1]+p[i]))
        end
    end
    Aitken(p)
end


#Modificar intervals per anar trobant a
A = [
[2.9,3.1],
[3.4,3.5],
[3.54,3.547],
[3.564,3.567],
[3.5685,3.569],
[3.56965,3.56970],
[3.56988,3.56994],
[3.56992,3.56995],
[3.56994,3.56999]
]
find_an(A,2)
ratios(A,2)