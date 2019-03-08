function value = randprb(A,P)

   x = sum(rand >= cumsum([0 P]));
   value = A(x);

end