#!/usr/bin/env python3
"""
Example Usage of Material Database
=================================

This script demonstrates how to use the material database to generate
different PyKO configurations with various material combinations.

Author: PyKO Enhanced
Date: 2024
"""

from material_database import (
    print_material_database, 
    get_available_materials,
    create_material_config_files
)

def main():
    """
    Demonstrate various material configurations.
    """
    print("🎯 PyKO Material Database Example Usage")
    print("=" * 50)
    
    # 1. Show available materials
    print("\n📋 Available Materials:")
    print("-" * 30)
    materials = get_available_materials()
    for i, material in enumerate(materials, 1):
        print(f"{i:2d}. {material}")
    
    # 2. Print detailed material database
    print("\n📊 Detailed Material Properties:")
    print("-" * 40)
    print_material_database()
    
    # 3. Example configurations
    print("\n🚀 Generating Example Configurations:")
    print("-" * 40)
    
    # Example 1: Aluminum on Copper (classic combination)
    print("\n1️⃣ Aluminum → Copper (Classic)")
    print("-" * 30)
    create_material_config_files(
        impactor_material='Aluminum',
        target_material='Copper',
        impactor_thickness=0.0001,  # 100 μm
        target_thickness=0.0002,    # 200 μm
        impactor_cells=100,
        target_cells=200,
        velocity_range=(100.0, 500.0),
        velocity_steps=5
    )
    
    # Example 2: Steel on Steel (same material)
    print("\n2️⃣ Steel → Steel (Same Material)")
    print("-" * 30)
    create_material_config_files(
        impactor_material='Steel',
        target_material='Steel',
        impactor_thickness=0.0001,  # 100 μm
        target_thickness=0.0002,    # 200 μm
        impactor_cells=100,
        target_cells=200,
        velocity_range=(200.0, 800.0),
        velocity_steps=4
    )
    
    # Example 3: Tungsten on Lead (high density contrast)
    print("\n3️⃣ Tungsten → Lead (High Density Contrast)")
    print("-" * 40)
    create_material_config_files(
        impactor_material='Tungsten',
        target_material='Lead',
        impactor_thickness=0.00005, # 50 μm (thinner due to high density)
        target_thickness=0.0003,    # 300 μm
        impactor_cells=50,
        target_cells=300,
        velocity_range=(500.0, 1500.0),
        velocity_steps=6
    )
    
    # Example 4: Titanium on Iron (aerospace application)
    print("\n4️⃣ Titanium → Iron (Aerospace)")
    print("-" * 30)
    create_material_config_files(
        impactor_material='Titanium',
        target_material='Iron',
        impactor_thickness=0.0001,  # 100 μm
        target_thickness=0.0002,    # 200 μm
        impactor_cells=100,
        target_cells=200,
        velocity_range=(300.0, 1000.0),
        velocity_steps=5
    )
    
    # Example 5: Gold on Aluminum (noble metal on light metal)
    print("\n5️⃣ Gold → Aluminum (Noble on Light)")
    print("-" * 35)
    create_material_config_files(
        impactor_material='Gold',
        target_material='Aluminum',
        impactor_thickness=0.00005, # 50 μm (thinner due to high density)
        target_thickness=0.0002,    # 200 μm
        impactor_cells=50,
        target_cells=200,
        velocity_range=(150.0, 600.0),
        velocity_steps=4
    )
    
    print("\n✅ All example configurations generated!")
    print("\n📁 Files created in: test18_a_multivelocity/")
    print("📝 Each configuration includes interface separation")
    print("\n🎯 Next Steps:")
    print("   1. Review the generated YAML files")
    print("   2. Modify parameters as needed")
    print("   3. Run PyKO simulations with: python pko-test18_a_multivelocity.py")
    print("   4. Analyze results and plots")

if __name__ == "__main__":
    main()
