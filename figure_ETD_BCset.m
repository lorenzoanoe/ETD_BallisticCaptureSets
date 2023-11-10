%% Script for generating figures of the ETD Ballistic Capture set

% Initialize input parameters:
% - system to be considered
%       sys=1 for Sun-Jupiter
%       sys=8 for Earth-Moon
% - grid stepsize
%       grid_scale = 400 for lower definition figures
%       grid_scale = 800 for higher definition figures -> datasets file size is large and figure generation is slow
%       grid_scale = 1400 for very high definition figures -> datasets are even larger. Not all the datasets for any GAMMA are provided
% - three-body energy parameter GAMMA
%       GAMMA = 0 for Asymptotic Ballistic Capture (ABC)
%       ...
%       until
%       GAMMA = GAMMA_max
%       GAMMA_max = 1.26 for Sun-Jupiter system
%       GAMMA_max = 1.36 for Earth-Moon system

%% Definition of the inputs
sys = 8 ;           % Earth-Moon system
grid_scale = 400 ;  % lower definition, faster generation
GAMMA = 0.84 ;

% Systems list
% 1 - Sun-Jupiter           mu = 9.54e-04
% 2 - Jupiter-Ganymede      mu = 7.80e-05   ABC in collision
% 3 - Jupiter-Europa        mu = 2.53e-05   ABC in collision
% 4 - Sun-Neptune           mu = 5.15e-05
% 5 - Sun-Mars              mu = 3.22e-07
% 6 - Sun-Earth             mu = 3.00e-06
% 7 - Sun-Psyche            mu = 1.15e-11
% 8 - Earth-Moon            mu = 1.22e-02
% 9 - Sun-Pluto             mu = 6.55e-09
% 10 - Pluto-Charon         mu = 1.04e-01   ABC in collision
% e.g. with sys = 8 the Earth-Moon system is chosen from the list

G = 6.674e-20 ;             % [km^3 / (kg*s^2)] universal gravitational constant
m_1vec(1) = 1.9885e+30 ;        m_2vec(1) = 1.89819e27 ;    R_refvec(1) = 778.4e6 ;       R_P_kmvec(1) = 71492 ;
m_1vec(2) = 1.89819E+27 ;       m_2vec(2) = 1.48E+23 ;      R_refvec(2) = 1070400 ;       R_P_kmvec(2) = 5262 ;
m_1vec(3) = 1.89819E+27 ;       m_2vec(3) = 4.80E+22 ;      R_refvec(3) = 670900 ;        R_P_kmvec(3) = 1561 ;
m_1vec(4) = 1.9885e+30 ;        m_2vec(4) = 1.0243E+26 ;    R_refvec(4) = 4498252900 ;    R_P_kmvec(4) = 49528 ;
m_1vec(5) = 1.9885e+30 ;        m_2vec(5) = 6.417E+23 ;     R_refvec(5) = 227900000 ;     R_P_kmvec(5) = 3402.4 ;
m_1vec(6) = 1.9885e+30 ;        m_2vec(6) = 398600/G ;      R_refvec(6) = 149.6e6 ;       R_P_kmvec(6) = 6378.0 ;
m_1vec(7) = 1.9885e+30 ;        m_2vec(7) = 2.287e19 ;      R_refvec(7) = 4.36921e11 ;    R_P_kmvec(7) = 278 ;
m_1vec(8) = 5.9724E+24 ;        m_2vec(8) = 7.348E+22 ;     R_refvec(8) = 384748 ;        R_P_kmvec(8) = 1737.4 ;
m_1vec(9) = 1.9885e+30 ;        m_2vec(9) = 1.303E+22 ;     R_refvec(9) = 5.87e9 ;        R_P_kmvec(9) = 2376.6 ;
m_1vec(10) = 1.303E+22 ;        m_2vec(10) = 1.52E+21 ;     R_refvec(10) = 19591.4 ;      R_P_kmvec(10) = 606 ;

mu_1 = m_1vec(sys) * G ;             % [km^3 / s^2] gravitational constant of m1
mu_2 = m_2vec(sys) * G ;             % [km^3 / s^2] gravitational constant of m2
R_ref = R_refvec(sys) ;
R_P = R_P_kmvec(sys) / R_ref ;       % [km] physical radius of m2

mu = mu_2 / (mu_1+mu_2) ;               % mass ratio of primaries
T = 2*pi*sqrt(R_ref^3/(mu_1+mu_2)) ;    % [s] orbital period
n = 1 / sqrt(R_ref^3/(mu_1+mu_2)) ;     % [rad/s] mean motion
r_Hill = (mu/3)^(1/3) ;                 % dimensionless radius of Hill's sphere

LP = lagrangePoints(mu) ; % position of the Lagrangian points in the baricentric synodic frame
x_CJL = LP(:,1) ;
y_CJL = LP(:,2) ;
CJ_L = (x_CJL).^2 + y_CJL.^2 + 2*(1-mu)./sqrt((x_CJL+mu).^2+y_CJL.^2) + 2*mu./sqrt((x_CJL-(1-mu)).^2+y_CJL.^2) ; % Jacobi constant on the Lagrangian points (velocity=0)

% Position of the Lagrangian points in the synodic frame centered in M1
L1 = LP(1) + mu ;
L2 = LP(2) + mu ;
L3 = LP(3) + mu ;
x_L4 = cosd(60) ;   % x_L4 = LP(4,1) + mu ;
y_L4 = sind(60) ;   % y_L4 = LP(4,2) ;
x_L5 = cosd(-60) ;
y_L5 = sind(-60) ;

% Creation of the grid
A = 3.5*r_Hill ;    % width of the grid generation
B = 4.5*r_Hill ;    % height of the grid generation
step_NA = round(r_Hill/grid_scale, 1, 'significant') ;
vett_x20 = [ 0 : -step_NA : -A, step_NA : step_NA : A*1.2 ] ;
vett_x20 = sort(vett_x20) ;
vett_y20 = [ 0 : -step_NA : -B, step_NA : step_NA : B ] ;
vett_y20 = sort(vett_y20) ;
l_x = length(vett_x20) ;
l_y = length(vett_y20) ;

% Upload Matlab structure for the capture set
if grid_scale ~= 400
    load( append("strsys", num2str(sys), "Gamma", num2str(GAMMA*100), "V16_step", num2str(grid_scale), ".mat") )
else
    load( append("strsys", num2str(sys), "Gamma", num2str(GAMMA*100), "V16.mat") )
end

% Delete BCs coming from collision (pre-BC collisions)
str.min_r2(str.eBW==-2) = [] ;
str.min_r2_1stRev(str.eBW==-2) = [] ;
str.x20(str.eBW==-2) = [] ;
str.y20(str.eBW==-2) = [] ;
str.sigma0(str.eBW==-2) = [] ;
str.NRev(str.eBW==-2) = [] ;
str.NRevALL(str.eBW==-2) = [] ;
str.TimeEn(str.eBW==-2) = [] ;
str.CasoEn(str.eBW==-2) = [] ;
str.Coll(str.eBW==-2) = [] ;
str.a(str.eBW==-2) = [] ;
str.e(str.eBW==-2) = [] ;
str.aBW(str.eBW==-2) = [] ;
str.eBW(str.eBW==-2) = [] ;

%% Figures
figure;
hold on

% Prograde BCs plot
if any(any(str.NRev>0))
    [~,idyPro]=ismembertol(str.y20(str.NRev>0), vett_y20) ;
    [~,idxPro]=ismembertol(str.x20(str.NRev>0), vett_x20) ;
    NRevPro = sparse(idyPro, idxPro, double(str.NRev(str.NRev>0)), l_y, l_x) ;

    contour_pro = contourc(vett_x20, vett_y20, full(NRevPro), [0.2,0.2]) ;
    len_integers = max( size( contour_pro(1,contour_pro(2,:) >= 1) ) ) ;
    contour_pro(:,contour_pro(2,:) >= 1 ) = zeros(2, len_integers).*[NaN; NaN] ;
    ProPoly = polyshape(transpose(contour_pro), 'KeepCollinearPoints', true) ;
    ProPoly = translate(ProPoly, 1, 0) ;
    
    h1 = plot(ProPoly, 'DisplayName','Prograde') ;
    set(h1,'linewidth', 1, 'EdgeColor', [0 0.4470 0.7410], 'FaceColor', [0 0.4470 0.7410], 'FaceAlpha', 0.5);
end

% Retrograde BCs plot
if any(any(str.NRev<0))
    [~,idyRetro]=ismembertol(str.y20(str.NRev<0), vett_y20) ;
    [~,idxRetro]=ismembertol(str.x20(str.NRev<0), vett_x20) ;
    NRevRetro = sparse(idyRetro, idxRetro, double(str.NRev(str.NRev<0)), l_y, l_x) ;

    contour_retro = contourc(vett_x20, vett_y20, full(abs(NRevRetro)), [0.2,0.2]) ;
    len_integers = max( size( contour_retro(1,contour_retro(2,:) >= 1) ) ) ;
    contour_retro(:,contour_retro(2,:) >= 1 ) = zeros(2, len_integers).*[NaN; NaN] ;
    RetroPoly = polyshape(transpose(contour_retro), 'KeepCollinearPoints', true) ; 
    RetroPoly = translate(RetroPoly, 1, 0) ;

    h2 = plot(RetroPoly, 'DisplayName','Retrograde') ;
    set(h2,'linewidth', 1, 'EdgeColor', [0.4660 0.6740 0.1880], 'FaceColor', [0.4660 0.6740 0.1880], 'FaceAlpha', 0.5);
end

% Early collisions plot
if any(any(str.Coll(str.NRev==0)>=0))
    str.Coll = str.Coll + 1 ;
    [~,idyColl]=ismembertol(str.y20(str.NRev==0), vett_y20) ;
    [~,idxColl]=ismembertol(str.x20(str.NRev==0), vett_x20) ;

    EarlyColl = sparse(idyColl, idxColl, double(str.Coll(str.NRev==0)), l_y, l_x) ;
    contour_Coll = contourc(vett_x20, vett_y20, full(EarlyColl), [0.2, 0.2]) ;
    len_integers = max( size( contour_Coll( 2,contour_Coll(2,:) >= 1 ) ) ) ;
    contour_Coll(:,contour_Coll(2,:) >= 1 ) = zeros(2, len_integers).*[NaN; NaN] ;
    EarlyCollPoly = polyshape(transpose(contour_Coll), 'KeepCollinearPoints', true) ;

    EarlyCollPoly = translate(EarlyCollPoly, [1 0]) ;
    h3 = plot(EarlyCollPoly, 'DisplayName','Early collisions') ;
    set(h3, 'FaceColor', 'k', 'EdgeColor', 'k'); % , 'FaceAlpha', 0.5);
end

% BC-collisions plot
index_CollBC = intersect(find(str.NRev~=0), find(str.Coll>0)) ;
if any(index_CollBC)
    [~,idyCollBC]=ismembertol(str.y20(index_CollBC), vett_y20) ;
    [~,idxCollBC]=ismembertol(str.x20(index_CollBC), vett_x20) ;

    NRevCollBC = sparse(idyCollBC, idxCollBC, double(str.NRev(index_CollBC)), l_y, l_x) ;
    contour_CollBC = contourc(vett_x20, vett_y20, full(abs(NRevCollBC)), [0.2, 0.2]) ;
    len_integers = max( size( contour_CollBC( 2,contour_CollBC(2,:) >= 1 ) ) ) ;
    contour_CollBC(:,contour_CollBC(2,:) >= 1 ) = zeros(2, len_integers).*[NaN; NaN] ;
    CollBC = polyshape(transpose(contour_CollBC), 'KeepCollinearPoints', true) ;
    CollBC = translate(CollBC, 1, 0) ;

    hCollBC = plot(CollBC, 'DisplayName', 'Post-BC collisions') ;
    set(hCollBC,'linewidth', 0.8, 'EdgeColor', [0.6 0.4470 0.7410], 'FaceColor', [0.6 0.4470 0.7410], 'FaceAlpha', 0.8);
end

% All collisions plot
% if any(any(str.Coll))
%     str.Coll = str.Coll + 1 ;
%     [~,idyColl]=ismembertol(str.y20(str.Coll>0), vett_y20) ;
%     [~,idxColl]=ismembertol(str.x20(str.Coll>0), vett_x20) ;
% 
%     Coll = sparse(idyColl, idxColl, double(str.Coll(str.Coll>0)), l_y, l_x) ;
%     contour_Coll = contourc(vett_x20, vett_y20, full(Coll), [0.2, 0.2]) ;
%     len_integers = max( size( contour_Coll( 2,contour_Coll(2,:) >= 1 ) ) ) ;
%     contour_Coll(:,contour_Coll(2,:) >= 1 ) = zeros(2, len_integers).*[NaN; NaN] ;
%     CollPoly = polyshape(transpose(contour_Coll), 'KeepCollinearPoints', true) ;
%     CollPoly = translate(CollPoly, 1, 0) ;
% 
%     h3 = plot(CollPoly, 'DisplayName','Collisions') ;
%     set(h3, 'FaceColor', 'k', 'EdgeColor', 'k'); % , 'FaceAlpha', 0.5);
% end

% Long-duration BCs plot
if any(any(abs(str.NRevALL)>=3))
    [~,idy3]=ismembertol(str.y20(abs(str.NRev)>=3), vett_y20) ;
    [~,idx3]=ismembertol(str.x20(abs(str.NRev)>=3), vett_x20) ;
    NRev3 = sparse(idy3, idx3, double(str.NRev(abs(str.NRev)>=3)), l_y, l_x) ;

    contour_3 = contourc(vett_x20, vett_y20, full(abs(NRev3)), [2.2,2.2]) ;
    len_integers = max( size( contour_3(1,contour_3(2,:) >= 1) ) ) ;
    contour_3(:,contour_3(2,:) >= 1 ) = zeros(2, len_integers).*[NaN; NaN] ;
    Poly3 = polyshape(transpose(contour_3), 'KeepCollinearPoints', true) ;
    Poly3 = translate(Poly3, 1, 0) ;
    
    hN3 = plot(Poly3, 'DisplayName','3+ Revolutions') ;
    set(hN3, 'linewidth', 1.5, 'EdgeColor', [0.8500 0.3250 0.0980], 'FaceColor', [0.8500 0.3250 0.0980], 'FaceAlpha', 0.8);
end

plot(1,0,'k+','markersize',10, 'DisplayName','Jupiter')
plot([L1,L2],[0,0],'ro', 'MarkerFaceColor', 'r','markersize',2, 'DisplayName','L1 and L2')
axis equal
axis square
axis tight
lgd = legend;
lgd.NumColumns = 2;
lgd.Location = "southoutside" ;
xlabel('X [LU]', 'Interpreter','latex')
ylabel('Y [LU]', 'Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','FontSize',16)
box on
set(gca,'FontSize',16)
set(gcf,'Position',[20 90 800 900])