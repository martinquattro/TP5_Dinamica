function A = getPoseRel(q, sub, supra)
    
    % Parametros del IRB140
    a =  [70 360 0 0 0 0];
    d = [0 0 0 380 0 0];
    alpha = [-90 0 90 -90 90 0] * pi/180;
    
    A = eye(4);
    if supra~=sub  % para los mismos indices se devuelve la identidad
        for i=(sub+1):supra
          A_aux = [ cos(q(i)), -sin(q(i)) * cos(alpha(i)), sin(q(i)) * sin(alpha(i)), a(i) * cos(q(i));
                    sin(q(i)), cos(q(i)) * cos(alpha(i)), -cos(q(i)) * sin(alpha(i)), a(i) * sin(q(i));
                    0 , sin(alpha(i)), cos(alpha(i)), d(i);
                    0 , 0 , 0, 1];

          A = A * A_aux;
        end
    end