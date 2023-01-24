using HCubature

function D_parallel(x,y) 
    a/32π^2(((a*L/2 - e^(-x*a/2)*sinh(y*a/2))*(a*L/2 + e^(x*a/2)*sinh(y*a/2)) - im*ε)^(-1) + 
            ((a*L/2 + e^(-x*a/2)*sinh(y*a/2))*(a*L/2 - e^(x*a/2)*sinh(y*a/2)) - im*ε)^(-1)) 
end

function F(x,y)
    exp(-(x^2 + y^2)/4σ^2 - im*x*Ω)*D_parallel(x,y)    
end