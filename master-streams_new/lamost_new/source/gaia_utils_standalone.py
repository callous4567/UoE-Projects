# Function to iterate over the provided list
from astroquery.gaia import Gaia
import pickle
username = 
password = 
Gaia.login(user=username,
           password=password)

def deltables(_remove_names):
    """
    Remove all these user tables with a try/except.
    :param _remove_names: list of names
    """
    for i in _remove_names:
        try:
            Gaia.delete_user_table(table_name=i, force_removal=True)
        except Exception as e:
            print(e)

def deljobs(_remove_jobs):
    """
    Removes *all* these job ids.
    :param _remove_jobs: array of job integers (unique to our username)
    """
    try:
        Gaia.remove_jobs(list(_remove_jobs))
    except Exception as e:
        print(e)

def upload_table(_table, _name):

    """
    :param _table: astropy table
    :param _name: unique name for table
    :return: table_id (identifier for the table = username._name)
    """

    Gaia.upload_table(upload_resource=_table, table_name=_name, table_description="For the thesis! Hallelujah!")
    table_id = "user_" + username + "." + _name

    return table_id

def obtain_dr3ids(_table_id, _match_table_name, dr2_column="dr2_source_id_local", magdif=0.1, angdist=1):

    """

    Will attempt to match DR3 source to the provided dr2_column in provided _table_id (remote.) Returns the table_id
    of the match, alongside deleting the job once complete/etc (will wait on job to avoid spamming external server.)

    **Deletion is hit-or-miss: recommend complete purge after all processes complete. **

    ***Relevant column for match to dr2 is "dr3_source_id"***

    :param _table_id: (stored by server) table id to crossmatch for DR3
    :param _match_table_name: The name (not id) for the job- use an rng or something.
    :param dr2_column: Column to use for DR2 ID. Default dr2_source_id_local
    :param magdif: maximum magnitude difference admissible for match. Default 0.1 mag
    :param angdist: maximum angular difference (mas) admissible for match. Default 1 mas
    :return: table_id of match
    """

    # Set up the job/query
    """
    dr2_source_id_local (equivalent to dr2_source_id)
    dr2_source_id (match found in DR2)
    dr3_source_id (the match to dr2_source_id from gaiadr3.dr2_neighbourhood match)
    """
    job = Gaia.launch_job_async(query=("""
                                          SELECT * 
                                          FROM {0}
                                          JOIN gaiadr3.dr2_neighbourhood 
                                              ON gaiadr3.dr2_neighbourhood.dr2_source_id = {0}.{1}
                                          WHERE ABS(magnitude_difference) < {2} 
                                          AND ABS(angular_distance) < {3}
                                          """).format(_table_id, dr2_column, magdif, angdist),
                                name=_match_table_name)
    Gaia.upload_table_from_job(job)
    Gaia.delete_user_table(table_name=_table_id, force_removal=True)
    Gaia.remove_jobs([job.jobid])
    job_table_id = "user_" + username + ".t" + str(job.jobid)

    return job_table_id

default_columns = ["dr2_source_id", "dr3_source_id", "ra", "dec",
                                              "parallax", "parallax_error", "pmra",
                                              "pmra_error", "pmdec", "pmdec_error", "radial_velocity",
                                              "radial_velocity_error"]
def obtain_download_data(_table_id, _source_join_name, _data_table_name, _savepath,
                         dr3_columns=default_columns, max_parallax=0.2, min_error_ratio=3):

    """
    Will take the provided table ID with dr3_source_id column, and obtain dr3_columns (list- see "default_columns" in
    gaia_utils_standalone, subject to max_parallax and min_error_ratio. Will download table and save it to savepath.
    Will return whether the job was successful (True/False bool.)

    :param _table_id: id of the table from obtain_dr3ids (has dr3 source ids)
    :param _source_join_name: misc name for job (get all dr3 columns)
    :param _data_table_name: misc name for job (get dr3_columns columns)
    :param _savepath: local save path for the obtained table
    :param dr3_columns: list [1,2,3...] of columns wanted, i.e. [pmra, pmdec, ...]
    :param max_parallax: maximum value of parallax to admit in mas. (1/parallax) in arcseconds = distance in parsecs
    :param min_error_ratio: minimum ratio of parallax to error: factor of three is minimum recommended.
    :return: bool on whether process was successful (True/False)
    """

    try:
        # Join our DR3_ID table with the main gaia_source table
        job = Gaia.launch_job_async(query=("""
                                              SELECT * 
                                              FROM {0}
                                              JOIN gaiadr3.gaia_source AS gaiadr3
                                                  ON gaiadr3.source_id = {0}.dr3_source_id
                                                  
                                              WHERE ABS(parallax) < {1}
                                              AND ABS(parallax_over_error) > {2} 
                                              """).format(_table_id, max_parallax, min_error_ratio),
                                    name=_source_join_name)
        Gaia.upload_table_from_job(job)
        Gaia.delete_user_table(table_name=_table_id)
        Gaia.remove_jobs([job.jobid])
        joined_id = "user_" + username + ".t" + job.jobid

        # Set up the string for selection col_1, col_2, col_3... ,col_n
        selection_string = dr3_columns[0]
        for i in range(1,len(dr3_columns)):
            selection_string += (",{}").format(dr3_columns[i])

        # Grab the columns we're interested in
        job = Gaia.launch_job_async(query=("""
                                              SELECT {0}
                                              FROM {1}
                                              """).format(selection_string, joined_id),
                                    name=_data_table_name)
        results = job.get_results()
        Gaia.remove_jobs([job.jobid])
        Gaia.delete_user_table(table_name=joined_id)

        with open(_savepath, "wb") as f:
            pickle.dump(file=f, obj=results)

        return True

    except Exception as e:

        with open(_savepath, "w") as f:
            f.write(str(e))

        return False

def up_obt_download_save(_table, _subtable_index, _final_savepath):

    """

    Upload _table using _subtable_index as a label to keep track of it, saving the final columns (under default_columns)
    to _final_savepath as an astropy table you can pickle_load. True/False return for if it was a success.

    :param _table: the subtable
    :param _subtable_index: the index
    :param _final_savepath: the savepath
    :return: True/False on whether it was a success or a fail.

    """

    # Generate name (using the provided subtable index for this subtable)
    uploaded_id = "local_upload_" + str(_subtable_index)
    id = upload_table(_table, uploaded_id)

    # Generate a random name for the dr3ids-matching job
    dr3id_match_name = "dr3match_" + uploaded_id
    id = obtain_dr3ids(id, dr3id_match_name)

    # Now, two random names for the jobs squirted out for the fullmatch, and the final columns we select
    fullmatch_name = uploaded_id + "_fullmatch"
    selected_columns_name = uploaded_id + "_columns"

    # Return whether job was success (True) or fail (False)
    return obtain_download_data(id, fullmatch_name, selected_columns_name, _savepath=_final_savepath)

# TODO:
"""
Upload table 
Crossmatch to Gaia DR3 within some distance threshold (i.e. ra_obs dec_obs columns)
Find best matches if possible within 2-3 arcseconds
Download table with appropriate columns we want
"""
def up_obt_download_save_TUPLEMASK(tuple):
    return up_obt_download_save(tuple[0], tuple[1], tuple[2])

