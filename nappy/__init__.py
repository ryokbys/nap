"""
`nappy` provides python scripts related to programs in the *nap* package.
"""
__all__ = ['pmd','napsys','atom','common','io','rdf','adf','msd','gaussian_smear',
           'util','vasp','mkcell','manipulate','units','database']

__author__ = 'RYO KOBAYASHI'
__version__ = '250522'

from . import napsys, io, pmd

_nappy_dir = '.nappy'

def get_nappy_dir():
    import os
    homedir = os.environ['HOME']
    return homedir + '/' +_nappy_dir


def get_git_hash():
    import subprocess
    import os
    try:
        # このスクリプトが存在するディレクトリを取得
        script_dir = os.path.dirname(os.path.abspath(__file__))
        
        # git -C オプションで、作業ディレクトリを指定してコマンド実行
        git_hash = subprocess.check_output(
            ['git', '-C', script_dir, 'rev-parse', '--short', 'HEAD'],
            stderr=subprocess.STDOUT # Git外で実行した際のエラー出力をキャプチャ
        ).decode('ascii').strip()
        
        return git_hash
    except (subprocess.CalledProcessError, FileNotFoundError):
        # Gitが入っていない環境や、.gitディレクトリが削除された場合
        return "unknown"

# print(f"Revision: {get_git_hash()}")
    
