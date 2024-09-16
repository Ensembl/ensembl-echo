import requests
import pymysql

desc = []

def get_filtered_subtree(tax_id, desc, max_retries=3, delay=5):
    url = f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/taxonomy/taxon/{tax_id}/filtered_subtree"
    try:
        response = requests.get(url)
        response.raise_for_status()  # Raises an HTTPError if the response status code is 4XX or 5XX
        data = response.json()
        #print(data)
        edges = data.get('edges')
        if edges:
            taxon_data = edges.get(str(tax_id))
            if taxon_data:
                visible_children = taxon_data.get('visible_children')
                #print(f"Visible children: {visible_children}")
                #print(len(visible_children))
                for i in visible_children:
                    #print(i)
                    child_data = edges.get(str(i))
                    if child_data:
                        children_status = child_data.get('children_status')
                        #print(f"Child Taxon ID: {i}, Children Status: {children_status}")
                        if children_status == 'NO_VISIBLE_CHILDREN':
                            desc.append(i)
                            #print("No visible status")
                        if children_status == 'HAS_MORE_CHILDREN':
                            print(i)
                            get_filtered_subtree(str(i), desc, max_retries=max_retries, delay=delay)
                    else:
                        print("no data")
            else:
                print("no taxon data")
        else:
            print("no edges")
    except requests.exceptions.RequestException as e:
        print(f"Attempt {attempt + 1} of {max_retries} failed: {e}")
        if attempt < max_retries - 1:
            print(f"Retrying in {delay} seconds...")
            time.sleep(delay)


all_desc = get_filtered_subtree('7100',desc, max_retries=3, delay=5)
print(f"All descendants: {len(desc)}")
