from flask import Flask, render_template, request, send_from_directory, jsonify
import json

app = Flask(__name__)


def style_sequence(sequence, query):
    if not query:
        return sequence

    lower_sequence = sequence.lower()
    lower_query = query.lower()

    styled_sequence = []
    start = 0

    while start < len(lower_sequence):
        match_index = lower_sequence.find(lower_query, start)
        if match_index == -1:
            styled_sequence.append(sequence[start:])
            break

        styled_sequence.append(sequence[start:match_index])

        styled_sequence.append(
            f"<span style='color: red;'>{sequence[match_index:match_index + len(query)]}</span>"
        )

        start = match_index + len(query)

    return "".join(styled_sequence)


@app.route("/", methods=["GET"])
def get_antibodies():
    antibodies_query = request.args.get("query", "").lower()
    json_file = "sequences.json"

    try:
        with open(json_file, "r") as file:
            all_antibodies = json.load(file)
    except FileNotFoundError:
        return "Error: sequences.json file not found.", 500

    for antibody in all_antibodies:
        antibody["styled_sequence"] = style_sequence(
            antibody["sequence"], antibodies_query
        )

    return render_template(
        "antibodies.html",
        antibodies_list=all_antibodies,
        antibodies_query=antibodies_query,
    )


@app.route("/compound-image", methods=["GET"])
def get_compound_image():
    smiles = request.args.get("smiles", "")
    try:
        id_value = smiles.split("=")[1]
        image_file = f"{id_value}.jpeg"
        return send_from_directory("static", image_file, mimetype="image/jpeg")
    except (IndexError, FileNotFoundError):
        return jsonify({"error": "Image not found"}), 404


if __name__ == "__main__":
    app.run(debug=True)
