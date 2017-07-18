#include <iostream>
#include <cmath> // For the pow() function
#include "io.h"
#include "/home/oliver/Programs/gems/gmml/includes/gmml.hpp"

using namespace MolecularModeling;

typedef std::vector<Atom*> AtomVector;
typedef std::vector<Residue*> ResidueVector;
typedef std::vector<std::string> StringVector;

double GetDistanceToAtom(Atom *A, Atom *otherAtom);
AtomVector ExtractProteinNonPolarHydrogens(AtomVector *A);
AtomVector ExtractGlycanNonPolarHydrogens(AtomVector *A);

int main(int argc, char *argv[])
{
    std::cout << "Hello World!" << std::endl;

    std::string working_Directory = Find_Program_Working_Directory();
    std::string installation_Directory = Find_Program_Installation_Directory();

    std::string parameterDirectory = installation_Directory + "/CurrentParams";

    //************************************************//
    // Details for loading in a PDB file              //
    //************************************************//

    std::vector<std::string> amino_libs, glycam_libs, other_libs, prep;
    amino_libs.push_back(parameterDirectory + "/amino12.lib");
    amino_libs.push_back(parameterDirectory + "/aminoct12.lib");
    amino_libs.push_back(parameterDirectory + "/aminont12.lib");

    glycam_libs.push_back(parameterDirectory + "/GLYCAM_amino_06j_12SB.lib");
    glycam_libs.push_back(parameterDirectory + "/GLYCAM_aminoct_06j_12SB.lib");
    glycam_libs.push_back(parameterDirectory + "/GLYCAM_aminont_06j_12SB.lib");

    other_libs.push_back(parameterDirectory + "/nucleic12.lib");
    other_libs.push_back(parameterDirectory + "/nucleic12.lib");
    other_libs.push_back(parameterDirectory + "/solvents.lib");

    prep.push_back(parameterDirectory + "/GLYCAM_06j-1.prep");

    std::string parameter_file_path = parameterDirectory + "/GLYCAM_06j.dat";
    //std::string ion_parameter_file_path = parameterDirectory + "/atomic_ions.lib";

    //************************************************//
    // Load PDB file                                  //
    //************************************************//

    Assembly receptor;
    receptor.BuildAssemblyFromPdbFile( (working_Directory + "/" + argv[1]), amino_libs, glycam_libs, other_libs, prep, parameter_file_path );
    receptor.BuildStructureByDistance();
    Assembly ligand;
    ligand.BuildAssemblyFromPdbFile( (working_Directory + "/" + argv[2]), amino_libs, glycam_libs, other_libs, prep, parameter_file_path );
    ligand.BuildStructureByDistance();

    // Get an AtomVector with just the hydogen atoms in assembly, saves looping time later.
    AtomVector temp_atoms = receptor.GetAllAtomsOfAssembly();
    AtomVector receptor_h_atoms = ExtractProteinNonPolarHydrogens(&temp_atoms);
    temp_atoms = ligand.GetAllAtomsOfAssembly();
    AtomVector ligand_h_atoms = ExtractGlycanNonPolarHydrogens(&temp_atoms);

    // Iterate through H atoms and calculate distance to atoms in target residue
    double distance = 0.0;
    double noe = 0.0, noe_total = 0.0, total_6VAH6 = 0.0, total_OME = 0.0, total_6VA = 0.0, total_0SA = 0.0, total_n4n6 = 0.0, total_n8n9 = 0.0;
    Atom *r_atom;
    Atom *l_atom;
    std::cout.precision(10);
    std::cout << std::fixed;
    for (AtomVector::iterator it = ligand_h_atoms.begin(); it != ligand_h_atoms.end(); ++it)
    {
        l_atom = *it;
        for (AtomVector::iterator itt = receptor_h_atoms.begin(); itt != receptor_h_atoms.end(); ++itt)
        {
            r_atom = *itt;
            distance = GetDistanceToAtom(r_atom, l_atom);
            //std::cout << "distance is " << distance << std::endl;
            if ( distance <= 6.0 && distance > 0.0 )
            {
                noe = (1 / pow (distance, 6) );
                noe_total += noe;
                std::cout << l_atom->GetResidue()->GetName() << ":" << l_atom->GetName() << " noe is " << noe << " with " << r_atom->GetResidue()->GetName() << ":" << r_atom->GetName() << std::endl;
            }
        }
        // Specific to my ligand, so I can add methyls together and other ambiguous peaks.
        if (l_atom->GetResidue()->GetName().compare("OME") == 0 )
        {
            total_OME += noe_total;
        }
        else if ( (l_atom->GetResidue()->GetName().compare("6VA") == 0 ) && (l_atom->GetName().find("M") != std::string::npos) )
        {
            total_6VA += noe_total;
        }
        else if ( (l_atom->GetResidue()->GetName().compare("6VA") == 0 ) && (l_atom->GetName().find("H6") != std::string::npos) )
        {
            total_6VAH6 += noe_total;
        }
        else if ( (l_atom->GetResidue()->GetName().compare("0SA") == 0 ) && (l_atom->GetName().find("M") != std::string::npos) )
        {
            total_0SA += noe_total;
        }
        else if ( (l_atom->GetResidue()->GetName().compare("0SA") == 0 ) && ( (l_atom->GetName().find("H8") != std::string::npos) || (l_atom->GetName().find("H9") != std::string::npos) ) )
        {
            total_n8n9 += noe_total;
        }
        else if ( (l_atom->GetResidue()->GetName().compare("0SA") == 0 ) && ( (l_atom->GetName().find("H4") != std::string::npos) || (l_atom->GetName().find("H6") != std::string::npos) ) )
        {
            total_n4n6 += noe_total;
        }
        else
        {
            std::cout << l_atom->GetResidue()->GetName() << ":" << l_atom->GetName() << " total noe is " << noe_total << std::endl;
        }
        noe_total = 0.0; // reset
    }

    std::cout << "OME:NAc total noe is " << total_OME << std::endl;
    std::cout << "6VA:NAc total noe is " << total_6VA << std::endl;
    std::cout << "6VA:H6 total noe is " << total_6VAH6 << std::endl;
    std::cout << "0SA:NAc total noe is " << total_0SA << std::endl;
    std::cout << "0SA:H4H6 total noe is " << total_n4n6 << std::endl;
    std::cout << "0SA:H8H9 total noe is " << total_n8n9 << std::endl;

    double tmp = (1 / pow (6, 6) );
    std::cout << "Smallest possible value is "<< tmp << std::endl;
    return 0;
}


double GetDistanceToAtom(Atom *A, Atom *otherAtom)
{
    double x = ( A->GetCoordinates().at(0)->GetX() - otherAtom->GetCoordinates().at(0)->GetX() );
    double y = ( A->GetCoordinates().at(0)->GetY() - otherAtom->GetCoordinates().at(0)->GetY() );
    double z = ( A->GetCoordinates().at(0)->GetZ() - otherAtom->GetCoordinates().at(0)->GetZ() );

    return sqrt( (x*x) + (y*y) + (z*z) );
}

AtomVector ExtractProteinNonPolarHydrogens(AtomVector *A)
{
    AtomVector nonPolarHydrogens;
    Atom *atom;
    for(AtomVector::iterator it = A->begin(); it != A->end(); ++it)
    {
        atom = *it;
        if (atom->GetName().find("H") != std::string::npos)
        {
            if (atom->GetName().find("HA") != std::string::npos) {nonPolarHydrogens.push_back(atom);}
            else if (atom->GetName().find("HB") != std::string::npos) {nonPolarHydrogens.push_back(atom);}
            else if (atom->GetName().find("HG2") != std::string::npos) {nonPolarHydrogens.push_back(atom);}
            else if ( (atom->GetName().find("HG") != std::string::npos) && (atom->GetResidue()->GetName().compare("SER") !=0) && (atom->GetResidue()->GetName().compare("THR") !=0)) {nonPolarHydrogens.push_back(atom);}
            else if ( (atom->GetName().compare("HD1") == 0) && (atom->GetResidue()->GetName().compare("HID") !=0) && (atom->GetResidue()->GetName().compare("HIS") !=0) ) {nonPolarHydrogens.push_back(atom);}
            else if (atom->GetName().compare("HD2") == 0) {nonPolarHydrogens.push_back(atom);}
            else if (atom->GetName().compare("HD3") == 0) {nonPolarHydrogens.push_back(atom);}
            else if ( (atom->GetName().compare("HE1") == 0) && (atom->GetResidue()->GetName().compare("TRP") !=0) ) {nonPolarHydrogens.push_back(atom);}
            else if ( (atom->GetName().compare("HE2") == 0) && (atom->GetResidue()->GetName().compare("HIE") !=0) && (atom->GetResidue()->GetName().compare("HIS") !=0) ) {nonPolarHydrogens.push_back(atom);}
            else if (atom->GetName().compare("HE3") == 0) {nonPolarHydrogens.push_back(atom);}
            else if ( (atom->GetName().find("HZ") != std::string::npos) && (atom->GetResidue()->GetName().compare("LYS") !=0) ) {nonPolarHydrogens.push_back(atom);}
            else if ( (atom->GetName().find("HD") != std::string::npos) && (atom->GetResidue()->GetName().compare("LEU") ==0) ) {nonPolarHydrogens.push_back(atom);}
            else if ( (atom->GetName().find("HD") != std::string::npos) && (atom->GetResidue()->GetName().compare("ILE") ==0) ) {nonPolarHydrogens.push_back(atom);}
            else if ( (atom->GetName().find("HE") != std::string::npos) && (atom->GetResidue()->GetName().compare("MET") ==0) ) {nonPolarHydrogens.push_back(atom);}
        }
    }
    return nonPolarHydrogens;
}

AtomVector ExtractGlycanNonPolarHydrogens(AtomVector *A)
{
    AtomVector nonPolarHydrogens;
    Atom *atom;
    for(AtomVector::iterator it = A->begin(); it != A->end(); ++it)
    {
        atom = *it;
        if (atom->GetName().find("H") != std::string::npos)
        {
            if ( (atom->GetName().find("O") != std::string::npos) || (atom->GetName().find("N") != std::string::npos) ) {}
            else if (atom->GetName().compare("CH3") == 0) {}
            else {nonPolarHydrogens.push_back(atom);}
        }
    }
    return nonPolarHydrogens;
}









