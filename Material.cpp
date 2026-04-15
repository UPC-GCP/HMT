// Imports
#include <iostream>
#include <vector> 
#include <json/json.h>

// Self-Imports
#include "Material.h"

Material::Material(Json::Value materials){
    
    // List
    vMat.resize(materials.size());

    // Store Materials
    for (Json::Value::ArrayIndex i = 0; i < materials.size(); i++){
        vMat[i].rho = materials[i]["rho"].asDouble();
        vMat[i].lambda = materials[i]["lambda"].asDouble();
        vMat[i].cp = materials[i]["cp"].asDouble();
        vMat[i].alpha = vMat[i].lambda / (vMat[i].rho * vMat[i].cp);
    }

}

void Material::setInitialConditions(double initTemp){

    // Initial Conditions
    T0 = initTemp;

}

