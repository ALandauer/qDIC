function [u_smooth] = equilb_smooth_qDIC(u,gridPoints,px2m,smoothing_params,plotting)
%Equilb smoothing for qDIC displacements

G = smoothing_params.G;
nu = smoothing_params.nu;

for ii = 1:2
    u{1}{ii} = px2m.*u{1}{ii};
    u{1}{ii} = inpaint_nans(u{1}{ii},4);
    gridPoints{ii} = px2m.*gridPoints{ii};
end


%----------- parameter definitions -----------

sizeU = size(u{1}{1});

ncoord = 2;
ndof = 2;
nnode = numel(u{1}{1});
nelem = (sizeU(1)-1)*(sizeU(2)-1);
maxnodes = 4;
nelnodes = 4*ones(1,nelem);

connect = image_connectivity(sizeU,nelem);

materialprops = [G, nu, 0]; %G, mu, plainstrain

% coords = [];
% for ii = 1:gridPoints{1}(end,1)
%     x_ = [1:gridPoints{2}(end,2);ii*ones(1,size(gridPoints{2}(:,2)))]';
%     coords = cat(1,coords,x_);
% end
% coords = coords'-1;

coords = [gridPoints{1}(:),gridPoints{2}(:)]';

% figure
% plotmesh(coords,ncoord,nnode,connect,nelem,nelnodes,'g');

%First, get the stiffness matrix needed for equilibrium

%find the full K matrix
K_ = globalstiffness(ncoord,ndof,nnode,coords,nelem,maxnodes,...
    nelnodes,connect,materialprops);

%trim out the boundary points where the equilibrium constraints are unknown

% m = 2*length(coords) - (2*sizeU(1)+2*sizeU(1));
% n = 2*length(coords);
% K = sparse([],[],[],m,n,0);
rows_kept = zeros(length(K_),1);
K_s = sparse(K_);
cnt = 0;

for ii = 1:numel(coords)
    
    %accept only the points that are not on the edges of u
    if coords(ii) ~= min(coords(1,:))
        if coords(ii) ~= min(coords(2,:))
            if coords(ii) ~= max(coords(1,:))
                if coords(ii) ~= max(coords(2,:))
                    cnt = cnt+1;
                    rows_kept(ii) = 1;
                end
            end
        end
    end
    
end
K = K_s(rows_kept==1,:);

% K = K_s(sizeU(1):end-sizeU(1),:);

%Second, solve the minimization

% figure,imagesc(K)
K_star = (K*K')\K;

% P = (speye(length(K_star)) - K'*K_star);
%rearrange u{}{} into a single vector
n = 0;
for jj = 1:size(u{1}{1},2)
    for kk = 1:size(u{1}{1},1)
        for ii = 1:2
            n = n+1;
            u_vect(n) = u{1}{ii}(kk,jj);
        end
    end
end

% P = (speye(length(K_star)) - K'*K_star);

%implement the eqn
u_smooth_vect = (speye(length(K_star)) - K'*K_star)*u_vect'; %eqn A21

%rearrange back into a matrix
n = 0;
for jj = 1:size(u{1}{1},2)
    for kk = 1:size(u{1}{1},1)
        for ii = 1:2
            n = n+1;
            u_smooth{1}{ii}(kk,jj) = u_smooth_vect(n)./px2m;
        end
    end
end

if plotting == 1
    figure
    subplot(1,2,1)
    surf(u{1}{1}./px2m)
    subplot(1,2,2)
    surf(u_smooth{1}{1})

    figure
    subplot(1,2,1)
    surf(u{1}{2}./px2m)
    subplot(1,2,2)
    surf(u_smooth{1}{2})
end

end

% ============== Req functions ===========================
%Functions are adapted from Prof. Bower's solidmechanics.org Matlab FE implementation

function Stif = globalstiffness(ncoord,ndof,nnode,coords,nelem,maxnodes,...
    nelnodes,connect,materialprops)
%
%   Assemble the global stiffness matrix
%

Stif = zeros(ndof*nnode,ndof*nnode);
lmncoord = zeros(ncoord,maxnodes);
%
%   Loop over all the elements
%
for lmn = 1:nelem
    %
    %   Extract coords of nodes, DOF for the current element
    %
    for a = 1:nelnodes(lmn)
        for i = 1:ncoord
            lmncoord(i,a) = coords(i,connect(a,lmn));
        end
    end
    n = nelnodes(lmn);
    kel = elstif(ncoord,ndof,n,lmncoord,materialprops);
    %
    %   Add the current element stiffness:the global stiffness
    %
    for a = 1:nelnodes(lmn)
        for i = 1:ndof
            for b = 1:nelnodes(lmn)
                for k = 1:ndof
                    rw = ndof*(connect(a,lmn)-1)+i;
                    cl = ndof*(connect(b,lmn)-1)+k;
                    Stif(rw,cl) = Stif(rw,cl) + kel(ndof*(a-1)+i,ndof*(b-1)+k);
                end
            end
        end
    end
end


end

function kel = elstif(ncoord,ndof,nelnodes,coord,materialprops)
%
%  Assemble the element stiffness
%
%    Arguments;
%
%      ncoord             No. coordinates (2 or 3 for 2D or 3D problem)
%      ndof               No. degrees of freedom per node (often ndof = ncoord)
%      nelnodes           No. nodes on the element
%      coords(i,a)        ith coord of ath node
%      materialprops      Material properties passed on:constitutive procedures

%   Local variables
%      npoints            No. integration points
%      xi(i,inpt)         ith local coord of integration point no. intpt
%      w(intpt)           weight for integration point no. intpt
%      N(a)               Shape function associated with ath node on element
%      dNdxi(a,i)         Derivative of ath shape function wrt ith local coord
%      dNdx(a,i)          Derivative of ath shape function wrt ith global coord
%      dxdxi(i,j)         Derivative of ith global coord wrt jth local coord
%      dxidx(i,j)         Derivative of ith local coord wrt jth global coord
%      det                Determinant of jacobian
%      strain(i,j)        strain_ij components
%      dsde(i,j,k,l)      Derivative of stress_ij with respect:strain_kl
%      kel(row,col)       Rows && cols of element stiffness
%
%
npoints = numberofintegrationpoints(ncoord,nelnodes);
dNdx = zeros(nelnodes,ncoord);
dxdxi = zeros(ncoord,ncoord);
kel = zeros(ndof*nelnodes,ndof*nelnodes);
%
%  Set up integration points && weights
%
xilist = integrationpoints(ncoord,nelnodes,npoints);
w = integrationweights(ncoord,nelnodes,npoints);
%
%  Loop over the integration points
%
for intpt = 1:npoints
    
    %     Compute shape functions && derivatives wrt local coords
    %
    for i = 1:ncoord
        xi(i) = xilist(i,intpt);
    end
    N = shapefunctions(nelnodes,ncoord,xi);
    dNdxi = shapefunctionderivs(nelnodes,ncoord,xi);
    
    
    %
    %     Compute the jacobian matrix && its determinant
    %
    for i = 1:ncoord
        for j = 1:ncoord
            dxdxi(i,j) = 0.;
            for a = 1:nelnodes
                dxdxi(i,j) = dxdxi(i,j) + coord(i,a)*dNdxi(a,j);
            end
        end
    end
    
    dxidx = inv(dxdxi);
    dt = det(dxdxi);
    %
    %     Convert shape function derivatives:derivatives wrt global coords
    %
    for a = 1:nelnodes
        for i = 1:ncoord
            dNdx(a,i) = 0.;
            for j = 1:ncoord
                dNdx(a,i) = dNdx(a,i) + dNdxi(a,j)*dxidx(j,i);
            end
        end
    end
    %
    %
    %     Compute the material tangent stiffness (d stress/d strain)
    %     ds/de is just C_ijkl for linear elasticity - this notation is used
    %     to allow extension to nonlinear problems
    %
    dsde = materialstiffness(ndof,ncoord,materialprops);
    %
    %     Compute the element stiffness
    %
    for a = 1:nelnodes
        for i = 1:ndof
            for b = 1:nelnodes
                for k = 1:ndof
                    row = ndof*(a-1)+i;
                    col = ndof*(b-1)+k;
                    for j = 1:ncoord
                        for l = 1:ncoord
                            kel(col,row) = kel(col,row) + ...
                                dsde(i,j,k,l)*dNdx(b,l)*dNdx(a,j)*w(intpt)*dt;
                        end
                    end
                end
            end
        end
    end
end

end

%
%====================== No. integration points =============================
%
%   Defines the number of integration points:be used for
%   each element type
%
function n = numberofintegrationpoints(ncoord,nelnodes)

if ncoord == 2
    if (nelnodes == 3)
        n = 1;
    end
    if (nelnodes == 6)
        n = 3;
    end
    if (nelnodes == 4)
        n = 4;
    end
    if (nelnodes == 8)
        n = 9;
    end
end

end
%
%====================== INTEGRATION POINTS ==================================
%
%   Defines positions of integration points
%
function xi = integrationpoints(ncoord,nelnodes,npoints)

xi = zeros(ncoord,npoints);
%
%    Rectangular element
%
if nelnodes(1) == 4
    if (npoints == 1)
        xi(1,1) = 0.;
        xi(2,1) = 0.;
    elseif (npoints == 4)
        xi(1,1) = -0.5773502692;
        xi(2,1) = xi(1,1);
        xi(1,2) = -xi(1,1);
        xi(2,2) = xi(1,1);
        xi(1,3) = xi(1,1);
        xi(2,3) = -xi(1,1);
        xi(1,4) = -xi(1,1);
        xi(2,4) = -xi(1,1);
    elseif (npoints == 9)
        xi(1,1) = -0.7745966692;
        xi(2,1) = xi(1,1);
        xi(1,2) = 0.0;
        xi(2,2) = xi(1,1);
        xi(1,3) = -xi(1,1);
        xi(2,3) = xi(1,1);
        xi(1,4) = xi(1,1);
        xi(2,4) = 0.0;
        xi(1,5) = 0.0;
        xi(2,5) = 0.0;
        xi(1,6) = -xi(1,1);
        xi(2,6) = 0.0;
        xi(1,7) = xi(1,1);
        xi(2,7) = -xi(1,1);
        xi(1,8) = 0.;
        xi(2,8) = -xi(1,1);
        xi(1,9) = -xi(1,1);
        xi(2,9) = -xi(1,1);
    end
end

end

%
%================= INTEGRATION WEIGHTS ==================================
%
%   Defines integration weights w_i
%
function w = integrationweights(ncoord,nelnodes,npoints)

w = zeros(npoints,1);

%
%    Rectangular element
%
if ncoord == 2
    if nelnodes(1) == 4
        if (npoints == 1)
            w(1) = 4.;
        elseif (npoints == 4)
            w = [1.,1.,1.,1.];
        elseif (npoints == 9 )
            w1D = [0.555555555,0.888888888,0.55555555555];
            for j = 1:3
                for i = 1:3
                    n = 3*(j-1)+i;
                    w(n) = w1D(i)*w1D(j);
                end
            end
        end
    end
end

end

% Calculate shape functions for the element
function N = shapefunctions(nelnodes,ncoord,xi)


N = zeros(nelnodes,1);
%
%    Rectangular element
%
if ncoord == 2
    N(1) = 0.25*(1.-xi(1))*(1.-xi(2));
    N(2) = 0.25*(1.+xi(1))*(1.-xi(2));
    N(3) = 0.25*(1.+xi(1))*(1.+xi(2));
    N(4) = 0.25*(1.-xi(1))*(1.+xi(2));
else
    disp('wrong number of coords')
end

end

%
%================= SHAPE FUNCTION DERIVATIVES ======================
%
function dNdxi = shapefunctionderivs(nelnodes,ncoord,xi)

dNdxi = zeros(nelnodes,ncoord);
%
%    Rectangular element
%
dNdxi(1,1) = -0.25*(1.-xi(2));
dNdxi(1,2) = -0.25*(1.-xi(1));
dNdxi(2,1) = 0.25*(1.-xi(2));
dNdxi(2,2) = -0.25*(1.+xi(1));
dNdxi(3,1) = 0.25*(1.+xi(2));
dNdxi(3,2) = 0.25*(1.+xi(1));
dNdxi(4,1) = -0.25*(1.+xi(2));
dNdxi(4,2) = 0.25*(1.-xi(1));

end
%
%================= Material Stiffness ==================================
%
%    Computes elasticity tensor C_{ijkl} = shear modulus and Poissons ratio
%    Currently coded either for plane strain, plane stress or general 3D.
%
function C = materialstiffness(ndof,ncoord,materialprops)

mu = materialprops(1);
nu = materialprops(2);

C = zeros(ndof,ncoord,ndof,ncoord);

%  planestrain = 0 => plane stress, planestrain = 1 => plane strain
planestrain = materialprops(3);

for i = 1:2
    for j = 1:2
        for k = 1:2
            for l = 1:2
                if (planestrain==1)
                    if (i==j && k==l)
                        C(i,j,k,l) = C(i,j,k,l)+2*mu*nu/(1-2*nu);
                    end
                else
                    if (i==j && k==l)
                        C(i,j,k,l) = C(i,j,k,l)+2*mu*nu/(1-nu);
                    end
                end
                if (i==l && k==j)
                    C(i,j,k,l) = C(i,j,k,l)+mu;
                end
                if (i==k && j==l)
                    C(i,j,k,l) = C(i,j,k,l)+mu;
                end
            end
        end
    end
end


end

function connect = image_connectivity(sizeU,nelem)
%For the simple case of a 2d image, define the element connectivity (4 nodes)

% 11- 12- 13- 14- 15
% | 5 | 6 | 7 | 8 |
% 6 - 7 - 8 - 9 - 10
% | 1 | 2 | 3 | 4 |
% 1 - 2 - 3 - 4 - 5


% % connect =
% %
% %      1     2     3     4     6     7     8     9
% %      2     3     4     5     7     8     9    10
% %      7     8     9    10    12    13    14    15
% %      6     7     8     9    11    12    13    14


connect_ = zeros(4,nelem);

cnt = 0;
for ii = 1:sizeU(2):nelem+1
    for jj = 1:sizeU(2)-1
        cnt = cnt+1;
        connect_(:,cnt) = [jj+ii-1,jj+ii,sizeU(2)+jj+ii,sizeU(2)+jj+ii-1];
    end
end

connect = connect_;

end


