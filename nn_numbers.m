%% Algorithm to identify hand written digits 0 - 10 from 20x20 pixel images. Use neural network with two hidden layers.

%Created by: Alex McDougal
%Last Editted: July 2018

%% Initialization
clear ; close all; clc

%% Setup the parameters
input_layer_size  = 400;  % 20x20 Input Images of Digits
hidden_layer_size = 25;   % 25 hidden units, layer 1
hidden_layer2_size=25;    % 25 hidden units, layer 2
num_labels = 10;          % 10 labels, from 1 to 10   

% Load Training Data
load('nn_data.mat');
m = size(X, 1);

%%Initialize neural network parameters

epsilon_init=0.12;
initial_Theta1=rand(hidden_layer_size,1+input_layer_size)*2*epsilon_init-epsilon_init;
initial_Theta2=rand(hidden_layer2_size,1+hidden_layer_size)*2*epsilon_init-epsilon_init;
initial_Theta3=rand(num_labels,1+hidden_layer2_size)*2*epsilon_init-epsilon_init;


% Unroll parameters
initial_nn_params = [initial_Theta1(:) ; initial_Theta2(:) ; initial_Theta3(:)];

%%Train Neural network

%Set parameters
options = optimset('MaxIter', 100);
lambda = 1;

% Abbreviate cost function
costFunction = @(p) nnCostFunction(p, ...
                                   input_layer_size, ...
                                   hidden_layer_size, hidden_layer2_size, ...
                                   num_labels, X, y, lambda);


[nn_params, cost] = fmincg(costFunction, initial_nn_params, options);

% Reshape theta 
Theta1 = reshape(nn_params(1:hidden_layer_size * (input_layer_size + 1)), ...
                 hidden_layer_size, (input_layer_size + 1));

Theta2 = reshape(nn_params((1 + (hidden_layer_size * (input_layer_size + 1))):(hidden_layer_size * (input_layer_size + 1))+(hidden_layer2_size*(hidden_layer_size+1))), ...
                 hidden_layer2_size, (hidden_layer_size + 1));

Theta3 = reshape(nn_params((1 + (hidden_layer_size * (input_layer_size + 1))+(hidden_layer2_size*(hidden_layer_size+1))):end), ...
                 num_labels, (hidden_layer2_size + 1));


		 
%%Predict training set accuracy
h1 = sigmoid([ones(m, 1) X] * Theta1');
h2 = sigmoid([ones(m, 1) h1] * Theta2');
h3 = sigmoid([ones(m, 1) h2] * Theta3');

[dummy, pred] = max(h3, [], 2);

fprintf('\nTraining Set Accuracy: %f\n', mean(double(pred == y)) * 100);


