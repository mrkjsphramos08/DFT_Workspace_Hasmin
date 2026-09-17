"""Measure in-plane layer distances in an LCO slab CIF.

Groups all species into horizontal z-layers, computes adjacent x-axis distances
within each layer and exports the combined table to ``lco_slab_distances.csv``.
"""
import pandas as pd

def analyze_full_lco_slab(cif_file_path, a_lattice=8.46351, output_csv="lco_slab_distances.csv", z_tol=0.03, y_tol=0.05):
    """
    Processes ALL species (Li, Co, O) in an LCO slab:
    1. Parses coordinates for all elements from the CIF file.
    2. Groups atoms into horizontal z-layers.
    3. Calculates adjacent x-axis distances (dx) for each layer.
    4. Automatically exports the combined dataset to CSV.
    """
    sites = []
    
    with open(cif_file_path, 'r') as f:
        lines = f.readlines()

    # 1. Parse atomic species and fractional coordinates
    for line in lines:
        parts = line.strip().split()
        if len(parts) >= 5 and not line.startswith('_') and not line.startswith('loop_') and not line.startswith('#'):
            # Detect species (Li, Co, O) and atom label
            species = ''.join([c for c in parts[0] if c.isalpha()]).capitalize()
            if not species or len(species) > 2:
                species = ''.join([c for c in parts[1] if c.isalpha()]).capitalize()
            
            label = parts[1] if len(parts) > 1 and parts[1] != parts[0] else parts[0]
            
            # Extract fractional x, y, z coordinates
            try:
                if len(parts) >= 6:
                    fx = float(parts[3].split('(')[0])
                    fy = float(parts[4].split('(')[0])
                    fz = float(parts[5].split('(')[0])
                    sites.append({'label': label, 'species': species, 'fx': fx, 'fy': fy, 'fz': fz})
            except ValueError:
                continue

    df = pd.DataFrame(sites).sort_values(by='fz').reset_index(drop=True)

    # 2. Group into distinct Z-layers across all species
    layers = []
    current_layer = [df.iloc[0]]
    for i in range(1, len(df)):
        if abs(df.iloc[i]['fz'] - current_layer[0]['fz']) <= z_tol:
            current_layer.append(df.iloc[i])
        else:
            layers.append(pd.DataFrame(current_layer))
            current_layer = [df.iloc[i]]
    if current_layer:
        layers.append(pd.DataFrame(current_layer))

    # 3. Calculate horizontal (a-axis) distances within each layer
    results = []
    for layer_idx, layer_df in enumerate(layers):
        avg_z = layer_df['fz'].mean()
        layer_df = layer_df.sort_values(by=['fy', 'fx']).reset_index(drop=True)
        
        for i in range(len(layer_df)):
            for j in range(i + 1, len(layer_df)):
                s1 = layer_df.iloc[i]
                s2 = layer_df.iloc[j]
                
                dfx = abs(s2['fx'] - s1['fx'])
                dfy = abs(s2['fy'] - s1['fy'])
                dfz = abs(s2['fz'] - s1['fz'])
                
                # Filter adjacent horizontal neighbors along x (fractional step ~ 0.33)
                if dfy <= y_tol and dfz <= z_tol and 0.20 <= dfx <= 0.40:
                    results.append({
                        'Layer_ID': f"Layer {layer_idx + 1} (z ≈ {avg_z:.2f})",
                        'Species_Type': f"{s1['species']} - {s2['species']}",
                        'Pair': f"{s1['label']} - {s2['label']}",
                        'Row_y': round((s1['fy'] + s2['fy']) / 2, 4),
                        'Δx_frac': round(dfx, 5),
                        'Distance_dx (Å)': round(dfx * a_lattice, 4)
                    })
                    
    results_df = pd.DataFrame(results)
    
    # 4. Save directly to output CSV
    results_df.to_csv(output_csv, index=False)
    print(f"Analysis complete! Processed {len(results_df)} total horizontal pairs.")
    print(f"Results saved to '{output_csv}'.")
    
    return results_df

df_all = analyze_full_lco_slab('LCO_slab_relaxed_3x2.cif', a_lattice=8.46351)