function [J grad] = nnCostFunction(nn_params, ...
                                   input_layer_size, ...
                                   hidden_layer_size, hidden_layer2_size, ...
                                   num_labels, ...
                                   X, y, lambda)
%NNCOSTFUNCTION Implements the cost function for a 3 layer neural network. Returns the cost function (J) and it's gradient (grad)
%Created by: Alex McDougal
%Last Editted: July 2018

% Reshape nn_params back into Theta1, Theta2, and Theta3
Theta1 = reshape(nn_params(1:hidden_layer_size * (input_layer_size + 1)), ...
                 hidden_layer_size, (input_layer_size + 1));

Theta2 = reshape(nn_params((1 + (hidden_layer_size * (input_layer_size + 1))):(hidden_layer_size * (input_layer_size + 1))+(hidden_layer2_size*(hidden_layer_size+1))), ...
                 hidden_layer2_size, (hidden_layer_size + 1));

Theta3 = reshape(nn_params((1 + (hidden_layer_size * (input_layer_size + 1))+(hidden_layer2_size*(hidden_layer_size+1))):end), ...
                 num_labels, (hidden_layer2_size + 1));


m = size(X, 1);
         
%Initialize
J = 0;
Theta1_grad = zeros(size(Theta1));
Theta2_grad = zeros(size(Theta2));
Theta3_grad = zeros(size(Theta3));



%calculate h
X_new=[ones(m,1) X];
z2=X_new*Theta1';
a2=sigmoid(z2);
a2_new=[ones(m,1) a2];
z3=a2_new*Theta2';
a3=sigmoid(z3);
a3_new=[ones(m,1) a3];
z4=a3_new*Theta3';;
h=sigmoid(z4);

%forward propagation
Delta1=0;
Delta2=0;
Delta3=0;
for i=1:m,
	%set up y vector
	y_i=zeros(num_labels,1);
	y_i(y(i))=1;

	h_i=h(i,:);
	%calc cost
	J=J+-1*log(h_i)*y_i - log(1-h_i)*(1-y_i);

	%calc delta (backwards propagation)
	d4=h_i'-y_i;
	d3=(Theta3(:,2:end)'*d4).*sigmoidGradient(z3(i,:))';

	d2=(Theta2(:,2:end)'*d3).*sigmoidGradient(z2(i,:))';
	
	Delta1=Delta1+d2*X_new(i,:);
	Delta2=Delta2+d3*a2_new(i,:);
	Delta3=Delta3+d4*a3_new(i,:);

endfor ;

J=J/m;
Theta1_grad=Delta1/m;
Theta2_grad=Delta2/m;
Theta3_grad=Delta3/m;

%add in regularization
reg1=Theta1.^2;
reg2=Theta2.^2;
reg3=Theta3.^2;
reg1(:,1)=0;
reg2(:,1)=0;
reg3(:,1)=0;
J=J+(lambda/(2*m))*(sum(reg1(:))+sum(reg2(:))+sum(reg3(:)));

reg4=Theta1;
reg5=Theta2;
reg6=Theta3;
reg4(:,1)=0;
reg5(:,1)=0;
reg6(:,1)=0;
Theta1_grad=Theta1_grad+(lambda/m)*reg4;
Theta2_grad=Theta2_grad+(lambda/m)*reg5;
Theta3_grad=Theta3_grad+(lambda/m)*reg6;



% Unroll gradients
grad = [Theta1_grad(:) ; Theta2_grad(:) ; Theta3_grad(:)];


end
