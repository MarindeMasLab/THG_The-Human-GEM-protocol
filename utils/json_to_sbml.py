from cobra.io import load_json_model, write_sbml_model

model_names = ["endoA_251209"]

for model_name in model_names:

    model = load_json_model(f"models/{model_name}.json")
    write_sbml_model(model, f"models/{model_name}.xml")