from read_my_file import *
def main():
    df = dataframe_process("report.pg_matrix 3.tsv")
    column_dict = column_group(df)
    result = dataframe_work(df, column_dict)

    print(result)

    return result

if __name__ == '__main__':
    result = main()
