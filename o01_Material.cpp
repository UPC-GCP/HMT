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

    // Main Initial Conditions
    Phi0 = initPhi;
    
}

void Material::setInitialConditions(double initPhi, Json::Value initVF){

    // Main Initial Conditions
    Phi0 = initPhi;
    
    // Velocity Field Initial Conditions
    VF0.resize(initVF.size());
    for (Json::Value::ArrayIndex i = 0; i < initVF.size(); i++){VF0[i] = initVF[i].asDouble();}

}
