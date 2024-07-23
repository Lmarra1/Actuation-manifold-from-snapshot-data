function Xd = isomap_decoder(G,X,Y,K)

% % Description of the input data
% X is the matrix with snapshots stored in columns
% Y is the matrix with points in the manifold stored in columns
% G is a matrix of points of the manifold to decode stored in columns
% K is the number of neighbours in the decoding algorithm based on KNN
% Xd are the decoded snapshots relative to G


%% Run the decoder
for i = size(G,2):-1:1
    
    % Read the i-th point in the manifold to be decoded
    yi = G(:,i);
    
    % Look for its K-nearest neighbours
    norm_yi = sqrt(sum((Y-yi).^2,1)); 
    [~,ind_neigb] = sort(norm_yi);
    
    % Store the points in latent and equivalent space of the first K
    % nearest neighbours
    ind_neigb = ind_neigb(1:K);
    neigb_y = Y(:,ind_neigb);
    neigb_x = X(:,ind_neigb);
    
    
    % Calculate matrix for function gradient estimation
    DY = (neigb_y(:,2:end) - neigb_y(:,1)).';
    DX = (neigb_x(:,2:end) - neigb_x(:,1)).';
    
    % Estimate the gradient of the decoding function
    f = pinv(DY)*DX;
    
    % Use first order Taylor expansion to obtain the point in the physical space
    Xd(:,i) = (neigb_x(:,1) + f.'*(yi-neigb_y(:,1)));
    
    
end



end