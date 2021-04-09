function [ S,I,R ] = SIR( S,I,R,N,x,y,c,m,v,e,r )

neighbor = zeros(x,y);
temp = 0;

%Calculate the sum needed by the model
if length(c) == 1
    for j = 2:y-1
        for i = 2:x-1

            if I(i-1,j-1) >= r
                temp = temp + (N(i-1,j-1)/N(i,j))*m*c*v*I(i-1,j-1);
            end
            if I(i-1,j) >= r
                temp = temp + (N(i-1,j)/N(i,j))*m*c*v*I(i-1,j);
            end
            if I(i-1,j+1) >= r
                temp = temp + (N(i-1,j+1)/N(i,j))*m*c*v*I(i-1,j+1);
            end
            if I(i,j-1) >= r
                temp = temp + (N(i,j-1)/N(i,j))*m*c*v*I(i,j-1);
            end
            if I(i,j+1) >= r
                temp = temp + (N(i,j+1)/N(i,j))*m*c*v*I(i,j+1);
            end
            if I(i+1,j-1) >= r
                temp = temp + (N(i+1,j-1)/N(i,j))*m*c*v*I(i+1,j-1);
            end
            if I(i+1,j) >= r
                temp = temp + (N(i+1,j)/N(i,j))*m*c*v*I(i+1,j);
            end
            if I(i+1,j+1) >= r
                temp = temp + (N(i+1,j+1)/N(i,j))*m*c*v*I(i+1,j+1);
            end
            if temp > 1
                temp = 1;
            end
            if temp < 0 
                temp = 0;
            end
            neighbor(i,j) = temp;
            temp = 0;
        end
    end
else 
    for j = 2:y-1
        for i = 2:x-1        
            if I(i-1,j-1) >= r
                temp = temp + (N(i-1,j-1)/N(i,j))*m*c(i-1,j-1)*v*I(i-1,j-1);
            end
            if I(i-1,j) >= r
                temp = temp + (N(i-1,j)/N(i,j))*m*c(i-1,j)*v*I(i-1,j);
            end
            if I(i-1,j+1) >= r
                temp = temp + (N(i-1,j+1)/N(i,j))*m*c(i-1,j+1)*v*I(i-1,j+1);
            end
            if I(i,j-1) >= r
                temp = temp + (N(i,j-1)/N(i,j))*m*c(i,j-1)*v*I(i,j-1);
            end
            if I(i,j+1) >= r
                temp = temp + (N(i,j+1)/N(i,j))*m*c(i,j+1)*v*I(i,j+1);
            end
            if I(i+1,j-1) >= r
                temp = temp + (N(i+1,j-1)/N(i,j))*m*c(i+1,j-1)*v*I(i+1,j-1);
            end
            if I(i+1,j) >= r
                temp = temp + (N(i+1,j)/N(i,j))*m*c(i+1,j)*v*I(i+1,j);
            end
            if I(i+1,j+1) >= r
                temp = temp + (N(i+1,j+1)/N(i,j))*m*c(i+1,j+1)*v*I(i+1,j+1);
            end
            if temp > 1
                temp = 1;
            end
            if temp < 0 
                temp = 0;
            end
            neighbor(i,j) = temp;
            temp = 0;
        end
    end
end

%Calculate the new S, I, R arrays
Smid = S - v*S.*I - S.*neighbor;
Imid = (1-e)*I + v*S.*I + S.*neighbor;
Rmid = R + e*I;

Smid(Smid < 0) = 0;
Imid(Imid > 1) = 1;
Imid(Imid <= 0.05) = 0;
Rmid(Rmid > 1) = 1;
Rmid(Rmid >= 0.95) = 1;

%Store the new arrays on S,I,R varibles
I = round(Imid,2);
S = round(Smid,2);
R = round(Rmid,2);
end

