#include<iostream>
#include<cmath>
using namespace std;
#define pi 3.1416
int main(){
    int ID;
    cout<<"Enter Your ID: ";
    cin>>ID;
    int R = ID % 100;
    double KVA = 25+(5*R);
    if (KVA <= 0){
        cout << "Error: KVA must be a positive value.\n";
        return 0;
    }
    double HV;
    double LV;
    double Temp = 75;
    // Check if R is even or odd and assign values
    if (R % 2 == 0) {
        HV = 6600;
        LV = 415;
    } else {
        HV = 11000;
        LV = 400;
    }
    int Number_of_legs = 3;
    double Voltage_per_turn = (sqrt((KVA*1000)/(Number_of_legs))) / (40);
    double CRGO_stell_lamination = 0.35;
    double Specific_magnetic_loading = 1.7;
    int Frequency = 50;
    double Net_cross_section_area_of_the_core = (Voltage_per_turn*1000000)/(4.44*Frequency*Specific_magnetic_loading);
    int Step = 7;
    double Iron_space_factor = 0.88;
    double Stacking_factor_for_lamination = 0.92;
    double Diameter_of_the_core = sqrt((4*Net_cross_section_area_of_the_core)/(Iron_space_factor*Stacking_factor_for_lamination*pi));
    Diameter_of_the_core = ((Diameter_of_the_core/10)+1)*10;
    Net_cross_section_area_of_the_core = Iron_space_factor*Stacking_factor_for_lamination*(pi/4)*Diameter_of_the_core*Diameter_of_the_core;
    Specific_magnetic_loading = (Voltage_per_turn*1000000)/(4.44*Frequency*Net_cross_section_area_of_the_core);

    cout << "\n=== CORE PARAMETERS ===\n";
    cout << "Voltage per turn: " << Voltage_per_turn << " V\n";
    cout << "Net cross-section area of core: " << Net_cross_section_area_of_the_core << " mm²\n";
    cout << "Diameter of core: " << Diameter_of_the_core << " mm\n";
    cout << "Specific magnetic loading: " << Specific_magnetic_loading << " T\n";

    double Window_space_factor = (10)/(30+(HV/1000));
    double Current_density = 2.5;
    double Window_area = (KVA*1000000000)/(3.33*Net_cross_section_area_of_the_core*Window_space_factor*Current_density*Specific_magnetic_loading*Frequency);

    cout << "\n=== WINDOW PARAMETERS ===\n";
    cout << "Window space factor: " << Window_space_factor << endl;
    cout << "Window area: " << Window_area << " mm²\n";

    double Voltage_per_phase_in_LV = LV/sqrt(3);
    int Number_of_turn_in_LV = Voltage_per_phase_in_LV / Voltage_per_turn;
    double Voltage_per_phase_in_HV = HV;
    int Number_of_turn_in_HV = (Voltage_per_phase_in_HV / Voltage_per_turn)*1.05;
    int Number_of_turn_in_HV_Normal = (Voltage_per_phase_in_HV / Voltage_per_turn);
    double Current_through_LV_winding = (KVA*1000)/ (sqrt(3)*LV);
    int Area_of_LV_conductor = Current_through_LV_winding/Current_density;

    cout << "\n=== WINDING PARAMETERS ===\n";
    cout << "Voltage per phase in LV: " << Voltage_per_phase_in_LV << " V\n";
    cout << "Number of turns in LV: " << Number_of_turn_in_LV << endl <<endl;
    cout << "Voltage per phase in HV: " << Voltage_per_phase_in_HV << " V\n";
    cout << "Number of turns in HV: " << Number_of_turn_in_HV << endl <<endl;
    cout << "Current through LV winding: " << Current_through_LV_winding << " A\n";
    cout << "Area of LV conductor: " << Area_of_LV_conductor << " mm²\n"<<endl;

    double Thinckness_of_LV_wire;
    double Width_of_LV_wire;
    // W and T constraints
    double W_min = 2, W_max = 16, T_min = 0.8, T_max = 5.6;
    const double tolerance = 0.01; // Allowable error for floating-point comparison

    bool found = false; // Flag to stop loops once a valid pair is found

    for (double W = W_min; W <= W_max; W += 0.1) { // Iterate over W with a step of 0.1
        for (double T = T_min; T <= T_max; T += 0.1) { // Iterate over T with a step of 0.1
            // Check if the area is approximately equal to 2 * T * W
            if (abs(2 * T * W - Area_of_LV_conductor) <= tolerance && (1.4 < W / T && W / T < 8)) {
                Width_of_LV_wire=(W * 10) / 10;
                Thinckness_of_LV_wire = (T * 10) / 10;
                found = true; // Set the flag to true
                break; // Exit the inner loop
            }
        }
        if (found) break; // Exit the outer loop if a pair is found
    }

    if (!found) {
        cout << "No valid values of W and T found for the given area.\n";
        return 0;
    }
    
    //Area_of_LV_conductor = 2*Thinckness_of_LV_wire*Width_of_LV_wire;

    cout<<"Thinckness of LV wire = "<<Thinckness_of_LV_wire<<endl;
    cout<<"Width of LV wire = "<<Width_of_LV_wire<<endl<<endl;

    double Current_through_HV_winding = (KVA*1000)/(3*HV);
    double Area_of_HV_conductor = Current_through_HV_winding / Current_density;
    double Diameter_of_the_conductor = sqrt((4*Area_of_HV_conductor)/(pi));
    double Copper_area_in_window = 2* ((Area_of_HV_conductor*Number_of_turn_in_HV)+(Area_of_LV_conductor*Number_of_turn_in_LV));
    Window_space_factor = (Copper_area_in_window/Window_area);

    cout<<"Current through HV winding = "<<Current_through_HV_winding<<endl;
    cout<<"Area of HV conductor = "<<Area_of_HV_conductor<<endl;
    cout<<"Diameter of the conductor = "<<Diameter_of_the_conductor<<endl<<endl;
    cout<<"Copper area in window = "<<Copper_area_in_window<<endl;
    cout<<"Window space factor = "<<Window_space_factor<<endl<<endl;

    int LV_layer = 2;

    double Paper_insulation_of_conductor = 0.25;
    double C_Thinckness_of_LV_wire = Thinckness_of_LV_wire + Paper_insulation_of_conductor;
    double C_Width_of_LV_wire = Width_of_LV_wire + Paper_insulation_of_conductor;
    double Turns_per_LV_layer = Number_of_turn_in_LV / LV_layer;
    double Height_of_LV_winding_in_window = Turns_per_LV_layer * C_Width_of_LV_wire;
    double Thickness_of_LV_coil = (2 * C_Thinckness_of_LV_wire) * LV_layer;
    double Distance_between_core_and_LV_coil = 3.5;
    double Mean_diameter_of_LV_coil = Diameter_of_the_core + 2 * Distance_between_core_and_LV_coil + (2 * C_Thinckness_of_LV_wire * LV_layer);
    double Mean_length_of_turn_of_LV_coil = 3.1416 * Mean_diameter_of_LV_coil;
    
    // Print LV Winding Values
    cout << "LV Winding Parameters:" << endl;
    cout << "Number of turns in LV winding: " << Number_of_turn_in_LV << endl;
    cout << "Corrected thickness of LV wire: " << C_Thinckness_of_LV_wire << " mm" << endl;
    cout << "Corrected width of LV wire: " << C_Width_of_LV_wire << " mm" << endl;
    cout << "Turns per LV layer: " << Turns_per_LV_layer << endl;
    cout << "Height of LV winding in window: " << Height_of_LV_winding_in_window << " mm" << endl;
    cout << "Thickness of LV coil: " << Thickness_of_LV_coil << " mm" << endl;
    cout << "Distance between core and LV coil: " << Distance_between_core_and_LV_coil << " mm" << endl;
    cout << "Mean diameter of LV coil: " << Mean_diameter_of_LV_coil << " mm" << endl;
    cout << "Mean length of turn of LV coil: " << Mean_length_of_turn_of_LV_coil << " mm" << endl; 

    // High Voltage (HV) Winding Parameters
    int HV_Layer = 11;

    int Splited_coils_in_HV_winding = 4;
    double Distance_between_LV_and_HV_coil = 12;
    double Changed_Diameter_of_the_conductor = Diameter_of_the_conductor + Paper_insulation_of_conductor;
    double Turns_in_each_coil = Number_of_turn_in_HV / Splited_coils_in_HV_winding;
    double Turns_per_HV_layer = (Turns_in_each_coil / HV_Layer);
    double Distance_in_between_every_HV_coil = 8;
    double Space_required_between_coils_and_core = 26;
    double Height_of_HV_coil_in_window = (Turns_per_HV_layer * Changed_Diameter_of_the_conductor * Splited_coils_in_HV_winding) +
                                         (Distance_in_between_every_HV_coil * (Splited_coils_in_HV_winding - 1));
    double Mean_diameter_of_HV_coil = Mean_diameter_of_LV_coil + Thickness_of_LV_coil +
                                      (2 * Distance_between_LV_and_HV_coil) +
                                      (Changed_Diameter_of_the_conductor * HV_Layer);
    double Mean_length_of_turn_of_HV_coil = 3.1416 * Mean_diameter_of_HV_coil;
    double Height_of_window_required = Height_of_HV_coil_in_window + (2 * Space_required_between_coils_and_core);
    double Window_width = Window_area / Height_of_window_required;

    // Print HV Winding Values
    cout << "\nHV Winding Parameters:" << endl;
    cout << "Number of turns in HV winding: " << Number_of_turn_in_HV << endl;
    cout << "Corrected diameter of HV conductor: " << Changed_Diameter_of_the_conductor << " mm" << endl;
    cout << "Turns in each coil: " << Turns_in_each_coil << endl;
    cout << "Turns per HV layer: " << Turns_per_HV_layer << endl;
    cout << "Height of HV coil in window: " << Height_of_HV_coil_in_window << " mm" << endl;
    cout << "Mean diameter of HV coil: " << Mean_diameter_of_HV_coil << " mm" << endl;
    cout << "Mean length of turn of HV coil: " << Mean_length_of_turn_of_HV_coil << " mm" << endl;
    cout << "Height of window required: " << Height_of_window_required << " mm" << endl;
    cout << "Window width: " << Window_width << " mm" << endl;

    // Percent Reactance
    double Vacuum_permeability = 0.00000125663;
    double Width_of_HV_winding = Changed_Diameter_of_the_conductor * HV_Layer;
    double Width_of_LV_winding = C_Thinckness_of_LV_wire * 2 * LV_layer;
    double Mean_height_of_the_coil = (Height_of_LV_winding_in_window + Height_of_HV_coil_in_window) / 2;
    double Lmt = (Mean_length_of_turn_of_HV_coil + Mean_length_of_turn_of_LV_coil) / 2;

    double Percent_reactance = ((2 * 3.1416 * Frequency * Vacuum_permeability * Lmt * Current_through_LV_winding *
                                Number_of_turn_in_LV * (Distance_between_LV_and_HV_coil +
                                ((Width_of_HV_winding + Width_of_LV_winding) / 3))) /
                                (Mean_height_of_the_coil * Voltage_per_turn * 1000)) * 100;

    // Print Percent Reactance
    cout << "\n## Percent_Reactance ##"<< endl;
    cout << "Lmt = "<<Lmt<<" mm"<<endl;
    cout << "Mean height of the coil = "<<Mean_height_of_the_coil<<" mm"<<endl;
    cout << "\nPercent Reactance: " << Percent_reactance << "%" << endl;

    // Percent Resistance
    double roh20 = 0.01724;
    double alpha20 = 0.00393;
    double roh = roh20 * (1 + (alpha20 * (Temp - 20)));
    double resistance_of_LV_winding = (roh * Mean_length_of_turn_of_LV_coil * Number_of_turn_in_LV) / (Area_of_LV_conductor * 1000);
    double resistance_of_HV_winding = (roh * Mean_length_of_turn_of_HV_coil * Number_of_turn_in_HV) / (Area_of_HV_conductor * 1000);
    double Transformation_ratio = HV / (LV / 1.7305);
    double Base_resistance = HV / Current_through_HV_winding;
    double Equivalent_resistance = resistance_of_HV_winding + (resistance_of_LV_winding * (Transformation_ratio * Transformation_ratio));
    double Percent_resistance = (Equivalent_resistance / Base_resistance) * 100;

    // Print Percent Resistance
    cout << "\n## Percent Resistance ##"<< endl;
    cout << "Resistance of LV winding per phase = "<< resistance_of_LV_winding << endl;
    cout << "Resistance of HV winding per phase = "<< resistance_of_HV_winding << endl;
    cout << "Ratio of Transformation = "<< Transformation_ratio << endl;
    cout << "Equivalent Resistance referred to HV winding per phase = "<< Equivalent_resistance << endl;
    cout << "\nPercent Resistance: " << Percent_resistance << "%" << endl;

    // Percent Impedance
    double Percent_impedance = sqrt((Percent_reactance * Percent_reactance) + (Percent_resistance * Percent_resistance));
    cout << "\n## Percent Impedance ##"<< endl;
    cout << "\nPercent Impedance: " << Percent_impedance << "%" << endl;\

    // Overall Dimensions of the Transformer
    double Width_of_largest_stamping = 0.95 * Diameter_of_the_core;
    double Distance_between_core_centers = Window_width + Diameter_of_the_core;
    double Height_of_yoke = Width_of_largest_stamping;
    double Depth_of_yoke = Width_of_largest_stamping;
    double Overall_height_of_the_frame = Height_of_window_required + 2 * Height_of_yoke;
    double Overall_width_of_the_frame = 2 * Distance_between_core_centers + Width_of_largest_stamping ;

    cout << "\n## Overall Dimensions of the Transformer ##"<< endl;
    cout << "Width of the window = "<< Window_width <<" mm"<< endl;
    cout << "Height of the window = "<< Height_of_window_required <<" mm"<< endl;
    cout << "Diameter of the circumscribing circle = "<< Diameter_of_the_core <<" mm"<<endl;
    cout << "Width of the largest stamping = "<< Width_of_largest_stamping <<" mm"<< endl;
    cout << "Distance between core centers = "<< Distance_between_core_centers <<" mm"<< endl;
    cout << "Height of yoke = "<< Height_of_yoke <<" mm"<< endl;
    cout << "Depth of yoke = "<< Depth_of_yoke <<" mm"<< endl;
    cout << "Overall height of frame = "<< Overall_height_of_the_frame <<" mm"<< endl;
    cout << "Overall width of frame = "<< Overall_width_of_the_frame <<" mm"<< endl;

    // Weight of Iron core and yoke assembly
    double Volume_of_the_core_and_yoke = Net_cross_section_area_of_the_core * ((Overall_width_of_the_frame * 2) + (Height_of_window_required * 3));
    double Weight_of_iron_per_kg = 7850;
    double Weight_of_core_and_yoke = (Volume_of_the_core_and_yoke * Weight_of_iron_per_kg) / (1000000000);

    cout << "\n## Weight of Iron core and yoke assembly ##"<<endl;
    cout << "Volume of the core and yoke = "<< Volume_of_the_core_and_yoke <<" mm3" << endl;
    cout << "Weight of iron = " << Weight_of_iron_per_kg <<" kg/m3"<< endl;
    cout << "Weight of core and yoke = " << Weight_of_core_and_yoke <<" kg"<< endl;

    // Core Loss
    double Core_loss_per_kg = 1.3;
    double Core_loss_in_Transformer = Weight_of_core_and_yoke * Core_loss_per_kg;

    cout << "\n## Core Loss in Transformer ##" << endl;
    cout << "Core loss per kg = " << Core_loss_per_kg << " W/kg" << endl;
    cout << "Core loss in transformer = " << Core_loss_in_Transformer << " W" << endl;

    // Calculation of Magnetizing Volt Amperes
    double Volt_ampere_per_kg = 10;
    double Magnetizing_VA_for_the_Transformer = Weight_of_core_and_yoke * Volt_ampere_per_kg;

    cout << "\n## Magnetizing Volt Amperes for Transformer ##" << endl;
    cout << "Volt ampere per kg = " << Volt_ampere_per_kg << " VA/kg" << endl;
    cout << "Magnetizing VA for the transformer = " << Magnetizing_VA_for_the_Transformer << " VA" << endl;

    // Weight of LV winding (per limb)
    double density_of_copper = 8.89;
    double Weight_of_LV_winding_per_limb = (density_of_copper * Number_of_turn_in_LV * Area_of_LV_conductor * Mean_length_of_turn_of_LV_coil) / (1000000);

    cout << "\n## Weight of LV Winding (per limb) ##" << endl;
    cout << "Density of copper = " << density_of_copper << " g/cm³" << endl;
    cout << "Weight of LV winding per limb = " << Weight_of_LV_winding_per_limb << " kg" << endl;

    // Weight of HV winding (per limb)

    double Weight_of_HV_winding_per_limb = (density_of_copper * Number_of_turn_in_HV * Area_of_HV_conductor * Mean_length_of_turn_of_HV_coil) / (1000000);
    double Weight_of_HV_winding_per_limb_Normal = (density_of_copper * Number_of_turn_in_HV_Normal * Area_of_HV_conductor * Mean_length_of_turn_of_HV_coil) / (1000000);

    cout << "\n## Weight of HV Winding (per limb) ##" << endl;
    cout << "Weight of HV winding per limb = " << Weight_of_HV_winding_per_limb << " kg" << endl;
    cout << "Weight of HV winding per limb (Normal) = " << Weight_of_HV_winding_per_limb_Normal << " kg" << endl;

    double Total_weight_of_copper_in_transformer_Normal = 3 * (Weight_of_LV_winding_per_limb + Weight_of_HV_winding_per_limb_Normal); 
    double Total_weight_of_copper_in_transformer = 3 * (Weight_of_LV_winding_per_limb + Weight_of_HV_winding_per_limb);

    cout << "\n## Total Weight of Copper in Transformer ##" << endl;
    cout << "Total weight of copper in transformer (Normal) = " << Total_weight_of_copper_in_transformer_Normal << " kg" << endl;
    cout << "Total weight of copper in transformer = " << Total_weight_of_copper_in_transformer << " kg" << endl;

    

    return 0;
}