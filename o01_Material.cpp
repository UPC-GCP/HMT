// Imports
#include <vector> 
#include <json/json.h>

// Self-Imports
#include "o01_Material.h"

Material::Material(Json::Value materials){
    
    // List
    vMat.resize(materials.size());

    // Store Materials
    for (Json::Value::ArrayIndex i = 0; i < materials.size(); i++){
        vMat[i].rho = materials[i]["rho"].asDouble();
		vMat[i].gamma = materials[i]["gamma"].asDouble();
        vMat[i].cp = materials[i]["cp"].asDouble();
    }

}

void Material::setInitialConditions(double initPhi){

    // Initial Conditions
    Phi0 = initPhi;

}
