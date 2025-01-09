from flask import (
    Flask,
    Response,
    json,
    redirect,
    render_template,
    request,
    send_from_directory,
)

app = Flask(__name__)


# @app.route("/", methods=["GET"])
# def root():
#     return render_template("base.html")


@app.route("/", methods=["GET"])
def get_antibodies():
    antibodies_query = request.args.get("query", "")
    all_antibodies = []
    with open("sequences.json", "r") as file:
        data = json.load(file)
        all_antibodies = data
    app.logger.info("all antibodies: ", all_antibodies)

    antibodies_matching_query = all_antibodies
    # if antibodies_query != "":
    #     match_indices_list = search_antibodies(antibodies_query, all_antibodies)
    #     antibodies_matching_query = [
    #         all_antibodies[i]
    #         for i, match_indices in enumerate(match_indices_list)
    #         if len(match_indices) > 0
    #     ]

    return render_template(
        "antibodies.html",
        antibodies_list=antibodies_matching_query,
        antibodies_query=antibodies_query,
    )


@app.route("/compound-image", methods=["GET"], endpoint="get_compound_image")
def get_compound_image():
    smiles = request.args.get("smiles", "example_antibody")
    image_file = "example_antibody.jpeg"
    # highlighted_substructure = request.args.get("highlight", "")
    # if highlighted_substructure == "":
    try:
        return send_from_directory("/static", image_file, mimetype="image/jpeg")
    except FileNotFoundError:
        return "Image not found", 404
    # else:
    #     return Response(
    #         get_highlighted_image(smiles, highlighted_substructure),
    #         mimetype="image/svg+xml",
    #     )
