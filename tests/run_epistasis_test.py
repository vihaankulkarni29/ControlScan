import pandas as pd
from src.controlscan.epistasis import MutationMatrix

def main():
    print("--- ControlScan Epistasis Engine Test ---\n")
    
    # Simulated Genomic Data: 6 E. coli isolates
    # Notice how acrA_T104A and acrB_N596H almost always appear together
    dummy_data = [
        {'sample_id': 'Strain_001', 'mutations': ['acrA_T104A', 'acrB_N596H', 'gyrA_S83L']},
        {'sample_id': 'Strain_002', 'mutations': ['acrA_T104A', 'acrB_N596H']},
        {'sample_id': 'Strain_003', 'mutations': ['acrA_T104A']}, # Missing acrB!
        {'sample_id': 'Strain_004', 'mutations': ['acrA_T104A', 'acrB_N596H', 'parC_S80I']},
        {'sample_id': 'Strain_005', 'mutations': ['gyrA_S83L', 'parC_S80I']}, # Completely different mechanism
        {'sample_id': 'Strain_006', 'mutations': ['acrA_T104A', 'acrB_N596H']},
    ]

    # Initialize the engine
    matrix_engine = MutationMatrix(dummy_data)

    # 1. Build the Matrix
    print("[1] Building Co-Occurrence Matrix...")
    co_matrix = matrix_engine.build_co_occurrence_matrix()
    print(co_matrix)
    print("\n")

    # 2. Calculate Exclusivity Probability
    print("[2] Calculating Epistatic Probabilities...")
    mut_A = 'acrA_T104A'
    mut_B = 'acrB_N596H'
    
    score = matrix_engine.calculate_exclusivity(mut_A, mut_B)
    print(f"Jaccard Co-occurrence Score for {mut_A} & {mut_B}: {score:.2f}")
    
    if score > 0.7:
        print(">> BIOLOGICAL CONCLUSION: High Epistasis Detected. These mutations are structurally co-dependent.")
    else:
        print(">> BIOLOGICAL CONCLUSION: Independent mutations.")

if __name__ == "__main__":
    main()
