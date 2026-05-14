function a = iscomplex(X)
a = ~(isreal(X))&&isnumeric(X);
end