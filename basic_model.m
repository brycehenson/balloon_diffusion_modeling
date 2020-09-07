

% this is a first go at simulating the deflation of a balloon as a function of time
% with multi gas species
% this is a usefull resource http://www.civilized.com/pdffiles/balloon.pdf
% but it only is used for a single gas

% the elasticity of balloons is rather complicated 
% https://www.jstor.org/stable/80395?seq=15#metadata_info_tab_contents
% https://www.cpp.edu/~tknguyen/che313/MeasurePermeability_of_Balloon_Material.docx
% https://www.balloonhq.com/faq/science.html
% http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.1000.5142&rep=rep1&type=pdf

%%

addpath('./lib/Core_BEC_Analysis/lib/') %add the path to set_up_project_path, this will change if Core_BEC_Analysis is included as a submodule
                  % in this case it should be 
set_up_project_path
hebec_constants

%%




syst_const=[];
syst_const.uninflated_radius=35e-3;
syst_const.uninflated_thickness=0.1e-3; %Pa 
syst_const.mod_elastic=0.00005;%0.1e9;
syst_const.outside_total_pressure=101e3; %Pa
% have a vector of the constant outside partial pressure
syst_const.outside_par_pressure=[0.78,0.21,0]*syst_const.outside_total_pressure; %Pa
syst_const.par_labels={'nitrogen','oxygen','helium'}; %Pa
syst_const.par_molecular_mass=[28.014,31.999,4.002]*1e-3; %kg/mol
syst_const.temperature=300;
syst_const.univ.atoms_per_mole=6.02214076e23;
syst_const.univ.boltzmann_constant=const.kb;
syst_const.univ.ideal_gas_const=8.31446261815324;


syst_const.diff_consts=20.4e-6*[0.16,0.45,0.65]; %from table 17 of  Permeability of Rubber to Gases Junius David Edwards, Samuel Fisher Pickering
% change units so that the equation on page 344 is in
% delta V= const* A t / d
% where V is in m^3, A is in m^2 t is in seconds and d is in m
% correct to time in seconds
syst_const.diff_consts=syst_const.diff_consts/60;
% correct to area in meters
syst_const.diff_consts=syst_const.diff_consts*(100^2);
% correct to thickness in meters
syst_const.diff_consts=syst_const.diff_consts/100;
% correct to volume in m^3
syst_const.diff_consts=syst_const.diff_consts*1e-6;
% convert delta V into delta moles 
syst_const.diff_consts=syst_const.diff_consts*syst_const.outside_total_pressure/(syst_const.temperature*syst_const.univ.ideal_gas_const);

% % convert delat moles into delta mass
% syst_const.diff_consts=syst_const.diff_consts*syst_const.par_molecular_mass;


% define the inital fill of the balloon by the number of moles of each species
% for a 27cm diameter balloon this is about 0.010 m^3
initial_tot_volume=0.010;
% saying that the inital pressure is about 2.6kPa (20mmhg) above atmosphere (https://www.quora.com/What-is-the-pressure-inside-a-normal-party-balloon)
% https://www.youtube.com/watch?v=PK1UK5-B5FQ
par_pressure=[0,0,syst_const.outside_total_pressure+2.6e3];
num_moles_inital=par_pressure*initial_tot_volume/(syst_const.temperature*syst_const.univ.ideal_gas_const);



%
syst_const.uninflated_sa=4*pi*(syst_const.uninflated_radius)^2;

%%
sum(num_moles_inital)

r_out=find_balloon_radius(sum(num_moles_inital),syst_const)

%%

function out=dndt_species(num_moles_now,syst_const)
    % first we must find the total pressure inside the balloon
    % first we find the p*v product with the ideal gas law
    tot_moles=sum(num_moles_now)
    
    pv_prod= tot_moles*syst_const.temperature*syst_const.univ.ideal_gas_const;
    
    %using the formula pv
    
    


end



function r_out=find_balloon_radius(total_moles,syst_const)

    pv_prod=total_moles*syst_const.temperature*syst_const.univ.ideal_gas_const;
    root_fun=@(r) r.^4 - (syst_const.uninflated_sa^(1/2)) * ...
                ( ...
                + 3 * pv_prod ./ (2*syst_const.mod_elastic*( (4*pi).^2 ) ) ...
                - syst_const.outside_total_pressure * r.^2 ./ (2*syst_const.mod_elastic*4*pi ) ...
                +1/((4*pi)^2)  ...
                );
%%
    stfig('root fun')
    rsamp=linspace(0,1,1e3)
    plot(rsamp,root_fun(rsamp))
%%
    r_out=fzero(root_fun,syst_const.uninflated_radius*10);
end

