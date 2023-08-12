def print_usage():
    # print("For now just guess")
    print("Usage:")
    msg = "python main.py --filter-soft|--filter-medium|--filter-hard <path_to_folder_with_vcf_files> <path_to_folder_to_save_results>"
    print(msg)

if __name__ == '__main__':
    print_usage()
