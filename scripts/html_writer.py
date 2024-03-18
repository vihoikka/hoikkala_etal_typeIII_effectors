import pandas as pd
import argparse
import os

# Parse arguments
parser = argparse.ArgumentParser(description='Generates an HTML table with hyperlinks')
parser.add_argument('-m', '--mastertable', help='mastertable', required=True)
parser.add_argument('-vf', '--viz_dir_full', help='viz_dir full', required=True)
parser.add_argument('-vr', '--viz_dir_relative', help='viz_dir rel', required=True)
parser.add_argument('-o', '--output_html', help='Output HTML file', required=True)
args = parser.parse_args()

# Parse arguments into variables
mastertable_path = args.mastertable
viz_dir_full = args.viz_dir_full
viz_dir_relative = args.viz_dir_relative
output_html = args.output_html

# Read .tsv file into dataframe
mastertable = pd.read_csv(mastertable_path, sep='\t')

# Add new columns with file paths to images and and CSV tables
mastertable['Visualisation_path'] = viz_dir_relative + "/" + mastertable['Locus'] + '/' + mastertable['Locus'] + '_viz.png'
mastertable['Visualisation_table_path'] = viz_dir_relative + "/" + mastertable['Locus'] + '/' + mastertable['Locus'] + '_merged.csv'

#tidy up the table by removing the following columns: Cas10, Cas5, Cas7, Cas10_GGDD, Cas10_GGDD_coord, Cas10_GGDD_seq, Cas10_GGDE_coord, Cas10_GGDE_seq, Cas10_HD, Cas10_HD_list, Cas10_DH, Cas10_HD_coord, Cas10_DH_coord, Cas10_coord, HD_E-value, unk_x, mem_x, locus_id, unk_y, unk01, sample
mastertable = mastertable.drop(columns=['Cas10', 
                                        'Cas5', 
                                        'Cas7', 
                                        'Cas10_GGDD', 
                                        'Cas10_GGDD_coord', 
                                        'Cas10_GGDD_seq', 
                                        'Cas10_GGDE_coord', 
                                        'Cas10_GGDE_seq', 
                                        'Cas10_HD', 
                                        'Cas10_HD_list', 
                                        'Cas10_DH', 
                                        'Cas10_HD_coord', 
                                        'Cas10_DH_coord', 
                                        'Cas10_coord', 
                                        'HD_E-value', 
                                        'unk_x', 
                                        'mem_x', 
                                        'locus_id', 
                                        'unk_y', 
                                        'ca3',
                                        'ca4',
                                        'ca6',
                                        'sam-amp',
                                        'NE_ca3',
                                        'NE_ca4',
                                        'NE_ca6',
                                        'NE_sam-amp'])

# also drop Cas10_GGDE, GGDD_E-value, GGDE_E-value, ca5, no_of_unknowns, unknown_proteins, NE_ca5, val, rng, ae1, crn1, crn2, crn3, csx15, csx16, csx20, has_ring_nuclease, ring_nuclease, fusion_components, fusion_protein
mastertable = mastertable.drop(columns=['Cas10_GGDE',
                                        'GGDD_E-value',
                                        'GGDE_E-value',
                                        'ca5',
                                        'no_of_unknowns',
                                        'unknown_proteins',
                                        'NE_ca5'])

#rename column con_ca3 to ca3 and so on
mastertable = mastertable.rename(columns={'con_ca3': 'ca3', 'con_ca4': 'ca4', 'con_ca6': 'ca6', 'con_sam_amp': 'sam-amp'})

mastertable_excel = mastertable.copy()

#remove columns Visualisation_path, Visualisation_table_path, Visualisation_link
mastertable_excel = mastertable_excel.drop(columns=['Visualisation_path', 'Visualisation_table_path'])

#save excel to same path as html but with .xlsx extension
excel_out = output_html.replace('.html', '.xlsx')
mastertable_excel.to_excel(excel_out, index=False)

# Generate a separate HTML page for each image showing the picture
for index, row in mastertable.iterrows():
    pic_path = os.path.join(viz_dir_relative, row['Locus'], row['Locus'] + '_viz.png')
    #create a html script to display the image instead of its path
    pic_html_script = f'<img src="{pic_path}">'
    mastertable.loc[index, "Visualisation_link"] = pic_html_script

# Replace path url to HTML hyperlinks
for column in ['Visualisation_path', 'Visualisation_table_path']:
     mastertable[column] = mastertable[column].apply(lambda x: f'<a href="{x}" target="_blank">Open</a>')

# Create HTML table
mastertable_to_html = mastertable.to_html(escape=False, index=False, classes = 'display dataTable" id = "table')

# Adding "toggle" buttons and "image" divs to the 'Visualisation_link' column 
#for index in mastertable.index:
#    img_url = mastertable.loc[index, "Visualisation_link"] # Considering Visualisation_link is a url to image
#    mastertable.loc[index, 'Visualisation_link'] = f'<img src="{img_url}">' #just show the image directly
    #mastertable.loc[index, 'Visualisation_link'] = f'<button id="btn_{index}" class="toggle-image" data-img="img_{index}">Toggle Image</button><div class="visual_link" id="img_{index}" style="display: none;"><img src="{img_url}"></div>'

# Generate the HTML file
with open(output_html, 'w') as file:
    file.write(f'''
<!DOCTYPE html>
<html>
<head>
        <title>CRISPR-Cas type III loci</title>
        <!-- include the CSS for DataTable -->
        <link href="https://cdn.datatables.net/v/dt/jq-3.7.0/jszip-3.10.1/dt-1.13.8/b-2.4.2/b-colvis-2.4.2/b-html5-2.4.2/cr-1.7.0/fh-3.4.0/rr-1.4.1/sb-1.6.0/datatables.min.css" rel="stylesheet">

        <!-- include the jQuery library -->
        <script
        src='https://code.jquery.com/jquery-3.7.1.min.js'
        integrity='sha256-/JqT3SQfawRcv/BIHPThkBvs0OEvtFFmqPF/lYI/Cxo='
        crossorigin="anonymous"></script>

        <!-- include the DataTable library -->
        <script src="https://cdn.datatables.net/v/dt/jq-3.7.0/jszip-3.10.1/dt-1.13.8/b-2.4.2/b-colvis-2.4.2/b-html5-2.4.2/cr-1.7.0/fh-3.4.0/rr-1.4.1/sb-1.6.0/datatables.min.js"></script>
        <script src="https://cdnjs.cloudflare.com/ajax/libs/pdfmake/0.2.7/pdfmake.min.js"></script>
        <script src="https://cdnjs.cloudflare.com/ajax/libs/pdfmake/0.2.7/vfs_fonts.js"></script>


        <style>
            tfoot input {{
            width: 100%;
            padding: 3px;
            box-sizing: border-box;
        }}
        </style>
               
        <style>
            .hide img {{display: none}}
        </style>

        <script>
        $(document).ready(function () {{
               
            // This is for the search boxes for each column
            $('#table thead tr')
                .clone(true)
                .addClass('filters')
                .appendTo('#table thead');
               
            var table = $('table').DataTable({{
                "paging":   true,  // Enable or disable table pagination.
                "ordering": true,  // Enable or disable ordering of columns.
                "info":     true,  // Enable or disable table information display field.
                "searching": true, // Enable or disable search function.
                "fixedHeader": true,
                "colReorder": true,
                "rowReorder": true,
                "dom": "Qlfrtip",
                "initComplete": function () {{
                var api = this.api();

                // For each column
                api
                .columns()
                .eq(0)
                .each(function (colIdx) {{
                    // Set the header cell to contain the input element
                    var cell = $('.filters th').eq(
                        $(api.column(colIdx).header()).index()
                    );
                    var title = $(cell).text();
                    $(cell).html('<input type="text" placeholder="' + title + '" />');

                    // On every keypress in this input
                    $(
                        'input',
                        $('.filters th').eq($(api.column(colIdx).header()).index())
                    )
                        .off('keyup change')
                        .on('change', function (e) {{
                            // Get the search value
                            $(this).attr('title', $(this).val());
                            var regexr = '({{search}})'; //$(this).parents('th').find('select').val();

                            var cursorPosition = this.selectionStart;
                            // Search the column for that value
                            api
                                .column(colIdx)
                                .search(
                                    this.value != ''
                                        ? regexr.replace('{{search}}', '(((' + this.value + ')))')
                                        : '',
                                    this.value != '',
                                    this.value == ''
                                )
                                .draw();
                        }})
                        .on('keyup', function (e) {{
                            e.stopPropagation();

                            $(this).trigger('change');
                            $(this)
                                .focus()[0]
                                .setSelectionRange(cursorPosition, cursorPosition);
                        }});
                }});
                }},
            }});
        }});
        </script>

</head>
<body>
    <h1>CRISPR-Cas type III loci</h1>
    <button id="button">SHOW/HIDE</button>
    <div id="tableContainer">
    {mastertable_to_html}
    </div>
    <script>
    var button = document.getElementById('button');
    var body = document.body;

    button.onclick = function() {{
        body.className = body.className == 'hide' ? '' : 'hide';
    }}
</script>
</body>
</html>
    ''')


mastertable_to_html = mastertable.to_html(escape=False, index=False, border=0)