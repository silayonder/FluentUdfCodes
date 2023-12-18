# FluentUdfCodes

The UDF codes of ads_source_c, ads_source_w, and ads_source_ht provide source terms for species transport, adsorption rate and energy equations for the computational domain. 
The macro used in diffusivityICMlong_f implements different diffusion coefficient values for different cell zones. 
Step boundary condition for species mass fraction is defined in inlet-stepBC_6000s. 
rhoYi_3Ddotproduct_fluid and rhoYi_3Ddotproduct_solid implement flux continuity for species mass fraction at the fluid-solid conjugate boundary based on the FVM scheme.
